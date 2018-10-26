#include "oneshot_dwdx.hpp"
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

void oneshot_dwdx(
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
	const int n_dvar = opt_opts.n_design_variables;
	const int n_elem = flo_opts.n_elem;
	const int n_resi = 3*n_elem;
	const int n_face = n_elem+1;

	if(n_elem!=x.size()-1) abort();
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
    for (int i = 1; i < n_elem; i++) {
        p = flo_opts.outlet_p;
        rho = p / (flo_opts.R * inlet_T);
        c = sqrt(gam * p / rho);
        u = c * flo_opts.inlet_mach;
        e = rho * (flo_opts.Cv * inlet_T + 0.5 * pow(u, 2));

        flow_data.W[i*3+0] = rho;
        flow_data.W[i*3+1] = rho * u;
        flow_data.W[i*3+2] = e;
    }
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


    std::vector<double> gradient_norm_list, timeVec;
    clock_t tic = clock();

    // Initialize B
    H.setIdentity();
    //H = H * 1.0;
    //if (opt_opts.exact_hessian != 0) {
    //    H = getAnalyticHessian(opt_opts.hessian_type, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, current_design);
    //    H = invertHessian(H);
    //}

	quasiOneD(x, area, flo_opts, &flow_data); // Start with converged flow
    double residual_norm = 0.0;
	for (int k = 0; k < n_elem; k++) {
		for (int i = 0; i < n_elem; i++) {
			residual_norm = residual_norm + pow(flow_data.residual[i*3+k], 2);
		}
	}
	residual_norm = sqrt(residual_norm);

    double current_cost = evalFitness(dx, flo_opts, flow_data.W, opt_opts);
    VectorXd gradient(n_dvar);
    VectorXd oldGrad(n_dvar); //BFGS
    gradient = getGradient(opt_opts.gradient_type, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, current_design);
    double gradient_norm = 0;
    for (int i = 0; i < n_dvar; i++)
        gradient_norm += pow(gradient[i], 2);
    gradient_norm = sqrt(gradient_norm);

    int iteration = 0;
    while(
		(iteration < opt_opts.opt_maxit) &&

		((gradient_norm > opt_opts.opt_tol) ||
		(residual_norm > flo_opts.flow_tol)) )
    {
		current_cost = evalFitness(dx, flo_opts, flow_data.W, opt_opts);
        iteration++ ;

        // Get flow update
		residual_norm = 1;
		for (int i = 0; i < 1 && residual_norm > 1e-2; i++) {
			stepInTime(flo_opts, area, dx, &flow_data);

			// Calculating the norm of the density residual
			residual_norm = 0;
			for (int k = 0; k < 1; k++) {
				for (int i = 0; i < n_elem; i++) {
					residual_norm = residual_norm + pow(flow_data.residual[i*3+k], 2);
				}
			}
			residual_norm = sqrt(residual_norm);
		}

		// Get adjoint update
		pIpW = evaldCostdW(opt_opts, flo_opts, flow_data.W, dx);
		pIpX = dCostdDes;

        //auto max_dt = max_element(std::begin(flow_data.dt), std::end(flow_data.dt));
		//SparseMatrix<double> dRdW = eval_dRdW_dRdX_adolc(flo_opts, area, flow_data);
		//pGpW = identity - (*max_dt)/dx[1] * dRdW;


		// Get design update
		//dAreadDes = evaldAreadDes(x, dx, current_design);
		//dCostdArea = evaldCostdArea(n_elem); //0
		//dCostdDes = dCostdArea.transpose() * dAreadDes; //0
		//pGpX = -(*max_dt)/dx[1] * evaldRdArea(flo_opts, flow_data) * dAreadDes;

		eval_dGdW_dGdX_adolc(x, flo_opts, flow_data, current_design, &pGpW, &pGpX);
		gradient = pIpX.transpose() + pIpW.transpose()*(pGpW*pGpX);
		//gradient = pIpX.transpose() + pIpW.transpose()*(pGpX);

		double step_size = 1e+0;

        if (opt_opts.descent_type == 1) {
            search_direction =  -gradient;
			design_change = step_size*search_direction;
        } else if (opt_opts.descent_type == 2) {
            if (iteration > 1) {
                H_BFGS = BFGS(H, oldGrad, gradient, design_change);
				double t = 1e-1;
                H = t*H_BFGS + (1-t)*H;
            }
            std::cout<<H<<std::endl;
            std::cout<<H.eigenvalues()<<std::endl;
            search_direction = -H * gradient;
			//design_change = step_size*search_direction;
			design_change = search_direction;
		}


		printf("%-15s %-15s %-15s\n", "Current Design", "Gradient","Search Direction");
		int step = 1;
		//if(n_dvar/8>=1) step = n_dvar/8;
		if(step!=1) printf("Only printing 1 out of %d variables\n", step);
		for (int i=0; i<n_dvar; i+=step) {
			printf("%15.5e %15.5e %15.5e\n", current_design.design_variables[i], gradient[i], design_change[i]);
		}

        //double initial_alpha = 1.0;
		//current_cost = linesearch_backtrack_unconstrained(
		//	initial_alpha, x, dx, search_direction, gradient, current_cost, flo_opts, opt_opts, &design_change, &flow_data, &current_design);

		for (int i=0; i<n_dvar; i++) {
			current_design.design_variables[i] += design_change[i];
		}
		area = evalS(current_design, x, dx);
        oldGrad = gradient;
		//gradient = getGradient(opt_opts.gradient_type, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, current_design);
		printf("Iteration: %d, Flow Residual %23.14e Cost Function %23.14e Gradient Norm: %23.14e \n",
			iteration, residual_norm, current_cost, gradient_norm);

        gradient_norm = 0;
        for (int i = 0; i < n_dvar; i++) {
            gradient_norm += pow(gradient[i], 2);
		}
        gradient_norm = sqrt(gradient_norm);
        gradient_norm_list.push_back(gradient_norm);

		clock_t toc = clock();
        double elapsed = (double)(toc-tic) / CLOCKS_PER_SEC;
        timeVec.push_back(elapsed);
        std::cout<<"Time: "<<elapsed<<std::endl;

		printf("Iteration: %d, Flow Residual %23.14e Cost Function %23.14e Gradient Norm: %23.14e \n", iteration, residual_norm, current_cost, gradient_norm);
        std::cout<<"End of Design Iteration: "<<iteration<<std::endl<<std::endl<<std::endl;
    }

    std::cout<<"Final Gradient:"<<std::endl;
    std::cout<<gradient<<std::endl;

    std::cout<<std::endl<<"Final Design:"<<std::endl;
    for (int i = 0; i < n_dvar; i++)
        std::cout<<current_design.design_variables[i]<<std::endl;


    double final_cost = evalFitness(dx, flo_opts, flow_data.W, opt_opts);
    std::cout<<"Fitness: "<<final_cost<<std::endl;

    outVec(constants.case_name, "OptConv.dat", "w", gradient_norm_list);
    outVec(constants.case_name, "OptTime.dat", "w", timeVec);


    return;
}

