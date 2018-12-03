#include "oneshot_adjoint.hpp"
#include "structures.hpp"
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<iomanip>
#include<Eigen/Dense>
#include<stdlib.h>//exit
#include"optimizer.hpp"
#include"quasiOneD.hpp"
#include"convert.hpp"
#include"timestep.hpp"
#include"fitness.hpp"
#include"grid.hpp"
#include"gradient.hpp"
#include"residuald1.hpp"
#include"residual_derivatives.hpp"
#include"fixed_point_derivatives.hpp"
#include"cost_derivative.hpp"
#include"parametrization.hpp"
#include"analyticHessian.hpp"
#include"output.hpp"
#include<time.h>
#include<stdlib.h>     /* srand, rand */

using namespace Eigen;
void oneshot_adjoint(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &dx,
	const struct Flow_options &flo_opts,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &initial_design)
{
	// **************************************************************************************************************************************
	// Initialize the flow
	const double gam = flo_opts.gam;
	const int n_elem = flo_opts.n_elem;
	const int n_resi = 3*n_elem;
	const int n_face = n_elem+1;
	if(n_elem!=(int)x.size()-1) abort();
	const int n_dvar = opt_opts.n_design_variables;
	class Flow_data<double> flow_data(n_elem);
    // Inlet flow properties
    double inlet_T = isenT(gam, flo_opts.inlet_total_T, flo_opts.inlet_mach);
    double p = isenP(gam, flo_opts.inlet_total_p, flo_opts.inlet_mach);
    double rho = p / (flo_opts.R * inlet_T);
    double c = sqrt(gam * p / rho);
    double u = flo_opts.inlet_mach * c;
    double e = rho * (flo_opts.Cv * inlet_T + 0.5 * pow(u, 2));

	flow_data.W[0*3+0] = rho;
	flow_data.W[0*3+1] = rho * u;
	flow_data.W[0*3+2] = e;
    // Flow properties initialization with outlet
    // State Vectors Initialization
    for (int i = 1; i < n_elem+2; i++) {
        p = flo_opts.outlet_p;
        rho = p / (flo_opts.R * inlet_T);
        c = sqrt(gam * p / rho);
        u = c * flo_opts.inlet_mach;
        e = rho * (flo_opts.Cv * inlet_T + 0.5 * pow(u, 2));

        flow_data.W[i*3+0] = rho;
        flow_data.W[i*3+1] = rho * u;
        flow_data.W[i*3+2] = e;
    }
	class Flow_data<double> flow_data_linesearch = flow_data;
	// **************************************************************************************************************************************
	// Initialize the adjoint
	VectorXd pIpW(n_resi);
    SparseMatrix<double> pGpW(n_resi, n_resi);
    SparseMatrix<double> pGpX(n_resi, n_dvar);
	const int n_stencil = 3;
	const int n_state = 3;
	pGpW.reserve(n_elem*n_state*n_state*n_stencil);

	VectorXd adjoint(n_resi);
	VectorXd adjoint_update(n_resi);
	adjoint.setZero();
	SparseMatrix<double> identity(n_resi, n_resi);
    identity.setIdentity();
	// **************************************************************************************************************************************
	// Initialize the design
	struct Design<double> current_design = initial_design;
	current_design.design_variables = initial_design.design_variables;
    std::vector<double> area = evalS(current_design, x, dx);

    MatrixXd dAreadDes(n_face, n_dvar);
	VectorXd dCostdArea(n_face);
	dCostdArea = evaldCostdArea(n_elem);
    VectorXd dCostdDes(n_dvar);
	dCostdDes.setZero();
    VectorXd pIpX(n_dvar);
	pIpX.setZero();


    MatrixXd H(n_dvar, n_dvar),
	         H_BFGS(n_dvar, n_dvar);
    VectorXd search_direction(n_dvar), design_change(n_dvar);


    double current_cost = evalFitness(dx, flo_opts, flow_data.W, opt_opts);

    VectorXd gradient(n_dvar);
    VectorXd oldGrad(n_dvar); //BFGS

    std::vector<double> gradient_norm_list, timeVec;
    clock_t tic = clock();

    // Initialize B
    H.setIdentity();
    H = H * 1.0e-0;
    //if (opt_opts.exact_hessian != 0) {
    //    H = getAnalyticHessian(x, dx, flow_data.W, area, designVar, hessian_type);
    //    H = invertHessian(H);
    //}

    double residual_norm = 1.0;
    double gradient_norm = 0;
    for (int i = 0; i < n_dvar; i++)
        gradient_norm += pow(gradient[i], 2);
    gradient_norm = sqrt(gradient_norm);
    int iteration = 0;

    // Design Loop
	quasiOneD(x, area, flo_opts, &flow_data); // Start with converged flow
	// Calculating the norm of the density residual

	//for (int i = 1; i < flo_opts.flow_maxit && residual_norm > flo_opts.flow_tol; i++) {
	//	stepInTime(flo_opts, area, dx, &flow_data);
	//	residual_norm = 0;
	//	for (int k = 0; k < n_elem; k++) {
	//		for (int i = 0; i < n_elem; i++) {
	//			residual_norm = residual_norm + pow(flow_data.residual[i*3+k], 2);
	//		}
	//	}
	//	residual_norm = sqrt(residual_norm);
	//}
	// Get adjoint update
	pIpW = evaldCostdW(opt_opts, flo_opts, flow_data.W, dx);
	auto max_dt = max_element(std::begin(flow_data.dt), std::end(flow_data.dt));


    SparseMatrix<double> dRdW = eval_dRdW_dRdX_adolc(flo_opts, area, flow_data);
	// G = w - dt/dx * R
	pGpW = identity - (*max_dt)/dx[1] * dRdW;


	eval_dGdW_dGdX_adolc(x, flo_opts, flow_data, current_design, &pGpW, &pGpX);

	double adjoint_update_norm = 1;
	for (int i = 0; i < 1 && adjoint_update_norm > 1e-12; i++) {
		adjoint = pIpW + pGpW.transpose()*adjoint;

		//const double step_size = 1.0;
		//adjoint_update = pIpW + pGpW.transpose()*adjoint - adjoint;
		//adjoint = (1-step_size)*adjoint - step_size*adjoint_update;
		//adjoint_update_norm = adjoint_update.norm();

		//VectorXd old_adj = adjoint;
		//adjoint = pIpW + pGpW.transpose()*adjoint;
		//adjoint_update = adjoint - old_adj;
		//adjoint_update_norm = adjoint_update.norm();
		VectorXd residual = adjoint - (adjoint.transpose()*pGpW).transpose() - pIpW;
		printf("Iteration: %d, Adjoint_update norm %23.14e\n", i, residual.norm());
	}
	VectorXd residual = adjoint - (adjoint.transpose()*pGpW).transpose() - pIpW;

	std::cout<< residual <<std::endl<<std::endl;

    //dRdW = eval_dRdW_dRdX_adolc(flo_opts, area, flow_data);
	//pGpW = identity - dRdW;
	//SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver0;
	//slusolver0.compute(identity.transpose() - pGpW.transpose());
	//if (slusolver0.info() != 0)
	//	std::cout<<"Factorization failed. Error: "<<slusolver0.info()<<std::endl;
	//VectorXd adjoint2 = slusolver0.solve(pIpW);

	//SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
	//slusolver1.compute(dRdW.transpose());
	//if (slusolver1.info() != 0)
	//	std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
	//VectorXd psi = slusolver1.solve(pIpW);
	//std::cout<<adjoint-adjoint2<<std::endl<<std::endl;
	//std::cout<<adjoint-psi<<std::endl<<std::endl;
	//std::cout<<psi-adjoint2<<std::endl<<std::endl;

    while(
		(iteration < opt_opts.opt_maxit) &&

		((gradient_norm > opt_opts.opt_tol) ||
		(residual_norm > flo_opts.flow_tol)) )
    {
		current_cost = evalFitness(dx, flo_opts, flow_data.W, opt_opts);
        iteration++ ;

        // Get flow update
		residual_norm = 1;
        flow_data.current_CFL = flo_opts.CFL_min;
		for (int i = 0; i < 10000 && residual_norm > 1e-12; i++) {
			stepInTime(flo_opts, area, dx, &flow_data);

			residual_norm = norm2(flow_data.residual);
		}

		// Get adjoint update
		pIpW = evaldCostdW(opt_opts, flo_opts, flow_data.W, dx);
        //auto max_dt = max_element(std::begin(flow_data.dt), std::end(flow_data.dt));
		//dRdW = eval_dRdW_dRdX_adolc(flo_opts, area, flow_data);
		//pGpW = identity - (*max_dt)/dx[1] * dRdW;

		eval_dGdW_dGdX_adolc(x, flo_opts, flow_data, current_design, &pGpW, &pGpX);

		double step_size = 1e-0;
		double adjoint_update_norm = 1;
		for (int i = 0; i < 10000 && adjoint_update_norm > 1e-9; i++) {

			//adjoint_update = pIpW + pGpW.transpose()*adjoint;
			//adjoint_update = step_size*(adjoint_update-adjoint);
			//adjoint_update_norm = (adjoint_update).norm();
			//adjoint = adjoint + adjoint_update;

			adjoint_update = pIpW + pGpW.transpose()*adjoint;
			adjoint_update_norm = (adjoint_update - adjoint).norm();
			adjoint = adjoint_update;
			//std::cout<<i<<" "<<adjoint_update_norm<<std::endl<<std::endl;
		}
		step_size = 1e+0;

		// Get design update
		dAreadDes = evaldAreadDes(x, dx, current_design);
		//dCostdArea = evaldCostdArea(n_elem); //0
		//dCostdDes = dCostdArea.transpose() * dAreadDes; //0
		pIpX = dCostdDes;
		//pGpX = -(*max_dt)/dx[1]*evaldRdArea(flo_opts, flow_data) * dAreadDes;
		gradient = pIpX + pGpX.transpose()*adjoint;

        if (opt_opts.descent_type == 1) {
            search_direction = -gradient;
			design_change = step_size*search_direction;
        } else if (opt_opts.descent_type == 2) {
            if (iteration > 1) {
                H_BFGS = BFGS(H, oldGrad, gradient, design_change);
				double t = 1e-3;
                H = t*H_BFGS + (1-t)*H;
                //H = H_BFGS;
				std::cout<<H.eigenvalues()<<std::endl;
            }
            search_direction = -H * gradient;
			design_change = search_direction;
		}


		printf("%-15s %-15s %-15s\n", "Current Design", "Gradient","Search Direction");
		int step = 1;
		if(step!=1) printf("Only printing 1 out of %d variables\n", step);
		for (int i=0; i<n_dvar; i+=step) {
			printf("%15.5e %15.5e %15.5e\n", current_design.design_variables[i], gradient[i], search_direction[i]);
		}

		bool do_linesearch = false;
		if (do_linesearch) {
			double c1 = 1e-4;
			double c_pk_grad = c1 * gradient.dot(search_direction);
			struct Design<double> new_design = current_design;
			std::vector<double> new_area = evalS(new_design, x, dx);
			flow_data_linesearch.W = flow_data.W;
			stepInTime(flo_opts, new_area, dx, &flow_data_linesearch);
			double new_cost = evalFitness(dx, flo_opts, flow_data_linesearch.W, opt_opts);
			while(new_cost > (current_cost + step_size * c_pk_grad))
			{
				step_size = step_size * 0.1;
				printf("Alpha Reduction: %e\n", step_size);
				if (step_size < 1e-10) {
					printf("Error. Can't find step size. Returning with tiny step.\n");
					break;
				}

				for (int i = 0; i < n_dvar; i++) {
					new_design.design_variables[i] = current_design.design_variables[i] + step_size * search_direction[i];
				}
				new_area = evalS(new_design, x, dx);
				flow_data_linesearch.W = flow_data.W;
				stepInTime(flo_opts, new_area, dx, &flow_data_linesearch);
				new_cost = evalFitness(dx, flo_opts, flow_data_linesearch.W, opt_opts);
				printf("new_cost: %e\n", new_cost);
				printf("current_cost + step_size/2.0 * c_pk_grad: %e\n", current_cost + step_size/ 2.0 * c_pk_grad);
			}
			for (int i = 0; i < n_dvar; i++) {
				//design_change[i] = new_design.design_variables[i] - current_design.design_variables[i];
				current_design.design_variables[i] = new_design.design_variables[i];
			}
		} else {
			for (int i=0; i<n_dvar; i++) {
				current_design.design_variables[i] += design_change[i];
			}
		}
		area = evalS(current_design, x, dx);
        oldGrad = gradient;

        gradient_norm = gradient.norm();
        gradient_norm_list.push_back(gradient_norm);


		clock_t toc = clock();
        double elapsed = (double)(toc-tic) / CLOCKS_PER_SEC;
        timeVec.push_back(elapsed);
        std::cout<<"Time: "<<elapsed<<std::endl;

		printf("Iteration: %d, Flow Residual %23.14e Cost Function %23.14e Gradient Norm: %23.14e Adjoint_update norm %23.14e\n",
			iteration, residual_norm, current_cost, gradient_norm, adjoint_update_norm);
        std::cout<<"End of Design Iteration: "<<iteration<<std::endl<<std::endl<<std::endl;
    }

    std::cout<<"Final Gradient:"<<std::endl;
    std::cout<<gradient<<std::endl;

    std::cout<<std::endl<<"Final Design:"<<std::endl;
    for (int i = 0; i < n_dvar; i++) {
        std::cout<<current_design.design_variables[i]<<std::endl;
	}
    outVec(constants.case_name, "final_design.dat", "w", gradient_norm_list);


    double final_cost = evalFitness(dx, flo_opts, flow_data.W, opt_opts);
    std::cout<<"Fitness: "<<final_cost<<std::endl;

    outVec(constants.case_name, "OptConv.dat", "w", gradient_norm_list);
    outVec(constants.case_name, "OptTime.dat", "w", timeVec);


    return;
}

