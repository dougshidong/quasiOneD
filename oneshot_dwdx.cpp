#include "structures.h"
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<iomanip>
#include<Eigen/Dense>
#include<stdlib.h>//exit
#include"optimizer.h"
#include"quasiOneD.h"
#include"timestep.h"
#include"fitness.h"
#include"grid.h"
#include"gradient.h"
#include"residuald1.h"
#include"cost_derivative.h"
#include"parametrization.h"
#include"analyticHessian.h"
#include"output.h"
#include<time.h>
#include<stdlib.h>     /* srand, rand */

using namespace Eigen;

void oneshot_dwdx(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &dx,
	const struct Flow_options &flo_opts,
	const struct Optimization_options &opt_opts,
	const struct Design &initial_design)
{
	// **************************************************************************************************************************************
	// Initialize the flow
	const double gam = flo_opts.gam;
	int n_elem = flo_opts.n_elem;
	if(n_elem!=x.size()) abort();
	int n_dvar = opt_opts.n_design_variables;
	struct Flow_data flow_data;
	flow_data.dt.resize(n_elem);
	flow_data.W.resize(3*n_elem);
	flow_data.W_stage.resize(3*n_elem);
	flow_data.fluxes.resize(3*(n_elem+1)); // At the faces
	flow_data.residual.resize(3*n_elem);
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
	VectorXd pIpW(3*n_elem);
	SparseMatrix<double> pGpW;
	VectorXd adjoint(3*n_elem);
	VectorXd adjoint_update(3*n_elem);
	adjoint.setZero();
	SparseMatrix<double> identity(3*n_elem, 3*n_elem);
    identity.setIdentity();
	// **************************************************************************************************************************************
	// Initialize the design
	struct Design current_design = initial_design;
	current_design.design_variables = initial_design.design_variables;
    std::vector<double> area = evalS(current_design, x, dx);

    MatrixXd dAreadDes(n_elem + 1, n_dvar);
	VectorXd dCostdArea(n_elem + 1);
	dCostdArea = evaldCostdArea(n_elem);
    VectorXd dCostdDes(n_dvar);
	dCostdDes.setZero();
    VectorXd pIpX(n_dvar);
	pIpX.setZero();
    MatrixXd pGpX(3*n_elem, n_dvar);


    MatrixXd H(n_dvar, n_dvar),
	         H_BFGS(n_dvar, n_dvar);
    VectorXd search_direction(n_dvar), design_change(n_dvar);


    double current_cost = evalFitness(dx, flo_opts, flow_data.W, opt_opts);

    VectorXd gradient(n_dvar);
    VectorXd oldGrad(n_dvar); //BFGS
    gradient = getGradient(opt_opts.gradient_type, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, current_design);

    std::vector<double> gradient_norm_list, timeVec;
    clock_t tic = clock();

    // Initialize B
    H.setIdentity();
    H = H * 1.0;
    if (opt_opts.exact_hessian != 0) {
        H = getAnalyticHessian(opt_opts.hessian_type, opt_opts.cost_function, x, dx, area, flo_opts, flow_data, opt_opts, current_design);
        H = invertHessian(H);
    }

    double residual_norm = 1.0;
    double gradient_norm = 0;
    for (int i = 0; i < n_dvar; i++)
        gradient_norm += pow(gradient[i], 2);
    gradient_norm = sqrt(gradient_norm);
    int iteration = 0;

    // Design Loop
	quasiOneD(x, area, flo_opts, &flow_data); // Start with converged flow
	// Calculating the norm of the density residual
	residual_norm = 0;
	for (int k = 0; k < n_elem; k++) {
		for (int i = 0; i < n_elem; i++) {
			residual_norm = residual_norm + pow(flow_data.residual[i*3+k], 2);
		}
	}
	residual_norm = sqrt(residual_norm);
    while(
		(iteration < opt_opts.opt_maxit) &&

		((gradient_norm > opt_opts.opt_tol) ||
		(residual_norm > flo_opts.flow_tol)) )
    {
		current_cost = evalFitness(dx, flo_opts, flow_data.W, opt_opts);
        iteration++ ;

        // Get flow update
		residual_norm = 1;
		for (int i = 1; i < 2 && residual_norm > 1e-2; i++) {
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
		pGpW = evaldRdW(area, flo_opts, flow_data);
        auto max_dt = max_element(std::begin(flow_data.dt), std::end(flow_data.dt));
		pGpW = identity - (*max_dt)*evaldRdW(area, flo_opts, flow_data) / dx[1];

		// Get design update
		dAreadDes = evaldAreadDes(x, dx, current_design);
		//dCostdArea = evaldCostdArea(n_elem); //0
		//dCostdDes = dCostdArea.transpose() * dAreadDes; //0
		pIpX = dCostdDes;
		pGpX = -evaldRdArea(flo_opts, flow_data) * dAreadDes;

		gradient = pIpX.transpose() + pIpW.transpose()*(pGpW*pGpX);
		//gradient = pIpX.transpose() + pIpW.transpose()*(pGpX);

		double step_size = 1e+2;

        if (opt_opts.descent_type == 1) {
            search_direction =  -gradient;
			design_change = step_size*search_direction;
        } else if (opt_opts.descent_type == 2) {
            if (iteration > 1) {
                H_BFGS = BFGS(H, oldGrad, gradient, design_change);
                H = H_BFGS;
            }
            std::cout<<H<<std::endl;
            search_direction = -H * gradient;
			design_change = step_size*search_direction;
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

