#include "structures.hpp"
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<iomanip>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<stdlib.h>//exit
#include"timestep.hpp"
#include"convert.hpp"
#include"solve_linear.hpp"
#include"optimizer.hpp"
#include"quasiOneD.hpp"
#include"residuald1.hpp"
#include"residuald2.hpp"
#include"output.hpp"
#include<time.h>
#include<stdlib.h>     /* srand, rand */

using namespace Eigen;

void second_order_flow(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &area,
	const struct Flow_options &flo_opts)
{
	const double gam = flo_opts.gam;
	int n_elem = flo_opts.n_elem;
	if(n_elem!=x.size()) abort();
    std::vector<double> dx(n_elem);
    for (int i = 1; i < n_elem-1; i++) {
        dx[i] =  (x[i] - x[i-1])/2  +  (x[i+1] - x[i])/2 ;
    }
    dx[0] = x[1] - x[0];
    dx[n_elem-1] = x[n_elem-1] - x[n_elem-2];
	// **************************************************************************************************************************************
	// Initialize the flow
    clock_t tic = clock();
    clock_t toc;
    double elapsed;
    std::vector <int> itV(flo_opts.flow_maxit/flo_opts.print_freq);
    std::vector<double> normV(flo_opts.flow_maxit/flo_opts.print_freq);
    std::vector<double> timeVec(flo_opts.flow_maxit/flo_opts.print_freq);
	class Flow_data<double> flow_data;
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
	VectorXd lagrange_mult(3*n_elem);
	lagrange_mult.setZero();

    // Evaluate R
    VectorXd dw(3*n_elem);
    // Evaluate R
    VectorXd old_rhs(3*n_elem);
    VectorXd rhs(3*n_elem);
    getDomainResi(flo_opts, area, &flow_data);
    // Evaluate dRdW
    SparseMatrix<double> dRdW = evaldRdW(area, flo_opts, flow_data);
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(area, flo_opts, flow_data);
    // Evaluate lambda*ddRdWdW
    SparseMatrix<double> Identity(3*n_elem, 3*n_elem);
    Identity.setIdentity();
    SparseMatrix<double> lambda_ddRdWdW(3*n_elem, 3*n_elem);
    lambda_ddRdWdW.setZero();
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        lambda_ddRdWdW += lagrange_mult(Ri) * ddRdWdW[Ri];
    }
    lambda_ddRdWdW.setIdentity();
    MatrixXd dense_lddRdWdW = MatrixXd(lambda_ddRdWdW);
    
	// **************************************************************************************************************************************
    // Initialize B
    MatrixXd H;
    H.setIdentity();
    H = H * 1.0;

    double residual_norm = 1.0;
    int iteration = 0;

	// Calculating the norm of the density residual
	residual_norm = 0;
	for (int k = 0; k < n_elem; k++) {
		for (int i = 0; i < n_elem; i++) {
			residual_norm = residual_norm + pow(flow_data.residual[i*3+k], 2);
		}
	}
	residual_norm = sqrt(residual_norm);
    double relax1 = 5e1;
    double relax2 = 5e1;
    double factor = 0.99;
    for (int i = 0; i < 2000; i++)
    {
		stepInTime(flo_opts, area, dx, &flow_data);
    }
    while(
		(iteration < flo_opts.flow_maxit) &&
		(residual_norm > flo_opts.flow_tol) )
    {
        iteration++ ;
        std::cout<<"Iteration "<<iteration <<"   NormR "<<std::setprecision(15)<<residual_norm<<std::endl;

        getDomainResi(flo_opts, area, &flow_data);
        for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
            rhs(Ri) = -flow_data.residual[Ri];
        }

        dRdW = evaldRdW(area, flo_opts, flow_data);
        dRdW += Identity*relax1;
        lagrange_mult = solve_linear(dRdW, rhs, 0, 1e-14);

        // Evaluate ddRdWdW
        ddRdWdW = evalddRdWdW(area, flo_opts, flow_data);
        // Evaluate lambda*ddRdWdW
        lambda_ddRdWdW.setZero();
        for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
            //ddRdWdW[Ri].setIdentity();
            lambda_ddRdWdW -= lagrange_mult(Ri) * ddRdWdW[Ri];
        }
        lambda_ddRdWdW += Identity*relax2;
        relax1 = relax1*factor;
        relax2 = relax2*factor;
        //dense_lddRdWdW = BFGS(dense_lddRdWdW, old_rhs, rhs, dw);
        //old_rhs = rhs;
        //VectorXcd eigval = dense_lddRdWdW.eigenvalues();
        //std::cout<<eigval<<std::endl;

        rhs.setZero();
        rhs = lagrange_mult.transpose()*dRdW;

        dw = solve_linear(lambda_ddRdWdW, rhs, 0, 1e-14);

        for (int i = 0; i < 3 * n_elem; i++) {
            flow_data.W[i] += dw(i);
        }

        // Calculating the norm of the density residual
        residual_norm = 0;
        for (int i = 0; i < n_elem; i++) {
            residual_norm = residual_norm + pow(flow_data.residual[i*3+0], 2);
		}
        residual_norm = sqrt(residual_norm);

        // Monitor Convergence
        if (iteration%flo_opts.print_freq == 0) {
            if (flo_opts.print_conv == 1) {
                std::cout<<"Iteration "<<iteration <<"   NormR "<<std::setprecision(15)<<residual_norm<<std::endl;
            }
            itV[iteration / flo_opts.print_freq - 1] = iteration;
            normV[iteration / flo_opts.print_freq - 1] = residual_norm;

            toc = clock();
            elapsed = (double)(toc-tic) / CLOCKS_PER_SEC;
            timeVec[iteration / flo_opts.print_freq - 1] = elapsed;

			if (flo_opts.print_solution == 1) {
				for (int k = 0; k < 3; k++) {
					std::cout<<"W"<<k + 1<<std::endl;
					for (int i = 0; i < n_elem; i++) {
						printf("%d %25.15e\n", i, flow_data.W[i*3+k]);
					}
				}
			}
        }
    }

	std::vector<double> p_vec(n_elem);
	for (int i = 0; i < n_elem; i++) {
		p_vec[i] = get_p(flo_opts.gam, flow_data.W[i*3+0], flow_data.W[i*3+1], flow_data.W[i*3+2]);
	}
    outVec(flo_opts.case_name, "current_pressure.dat", "w", p_vec);
    std::cout<<"Flow iteration = "<<iteration<<"   Density Residual = "<<residual_norm<<std::endl;
}
