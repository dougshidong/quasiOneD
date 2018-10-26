// Time Stepping Schemes
#include "timestep.hpp"
#include "structures.hpp"
#include<vector>
#include<math.h>
#include<iostream>
#include<Eigen/Core>
#include<Eigen/LU>
#include<Eigen/Sparse>
#include "flux.hpp"
#include "convert.hpp"
//#include "residuald1.hpp"
#include "boundary_conditions.hpp"
#include<adolc/adolc.h>
#include"adolc_eigen.hpp"

#include <cmath>
#include <complex>

// Domain Residual R = FS_i+1/2 - FS_i-1/2 - Qi
template<typename dreal>
void getDomainResi( 
	const struct Flow_options &flo_opts,
	const std::vector<dreal> &area,
    class Flow_data<dreal>* const flow_data)
{
	inletBC(flo_opts, flow_data);
	outletBC(flo_opts, flow_data);
	getFlux(flo_opts, flow_data->W, &(flow_data->fluxes));
	const int n_elem = flo_opts.n_elem;
    for (int i_state = 0; i_state < 3; i_state++) {
        for (int i_res = 1; i_res < n_elem+1; i_res++) {
            const int kim = (i_res - 1) * 3 + i_state;
            const int ki  = (i_res + 0) * 3 + i_state;
            flow_data->residual[ki] = flow_data->fluxes[ki] * area[i_res] - flow_data->fluxes[kim] * area[i_res-1];

			if (i_state==1) { // Source Term
				dreal p = get_p(flo_opts.gam, flow_data->W[i_res*3+0], flow_data->W[i_res*3+1], flow_data->W[i_res*3+2]);
				flow_data->residual[ki] -= p * (area[i_res] - area[i_res-1]);
			}
        }
    }
}
template void getDomainResi<double>( const struct Flow_options &flo_opts, const std::vector<double> &area, class Flow_data<double>* const flow_data);
template void getDomainResi<adouble>( const struct Flow_options &flo_opts, const std::vector<adouble> &area, class Flow_data<adouble>* const flow_data);
template void getDomainResi<std::complex<double>>( const struct Flow_options &flo_opts, const std::vector<std::complex<double>> &area, class Flow_data<std::complex<double>>* const flow_data);

template<typename dreal>
void evaluate_dt( 
	const struct Flow_options &flo_opts,
    const std::vector<double> &dx,
    class Flow_data<dreal>* const flow_data) {

	dreal new_CFL = flow_data->current_CFL * pow((flow_data->old_residual_norm / flow_data->current_residual_norm), flo_opts.CFL_ramp);
	new_CFL = fmax(new_CFL, flo_opts.CFL_min);
	new_CFL = fmin(new_CFL, flo_opts.CFL_max);
	flow_data->current_CFL = new_CFL;
	//std::cout<<flow_data->current_CFL<<std::endl;

	// Calculate Time Step
	int n_elem = flo_opts.n_elem;
	const int first_cell = 0;
	const int last_cell = n_elem+1;
	for (int i = 1; i < n_elem+1; i++) {
		const dreal u = flow_data->W[i*3+1] / flow_data->W[i*3+0];
		const dreal c = get_c(flo_opts.gam, flow_data->W[i*3+0], flow_data->W[i*3+1], flow_data->W[i*3+2]);
		flow_data->dt[i] = (flow_data->current_CFL * dx[i]) / fabs(u + c);
	}
    flow_data->dt[first_cell] = flow_data->dt[first_cell+1];
    flow_data->dt[last_cell] = flow_data->dt[last_cell-1];
}
template void evaluate_dt( const struct Flow_options &flo_opts, const std::vector<double> &dx, class Flow_data<double>* const flow_data);
template void evaluate_dt( const struct Flow_options &flo_opts, const std::vector<double> &dx, class Flow_data<adouble>* const flow_data);

template<typename dreal>
void EulerExplicitStep(
	const struct Flow_options &flo_opts,
    const std::vector<dreal> &area,
    const std::vector<double> &dx,
    class Flow_data<dreal>* const flow_data);

template<typename dreal>
void lusgs(
	const struct Flow_options &flo_opts,
    const std::vector<dreal> &area,
    const std::vector<double> &dx,
    class Flow_data<dreal>* const flow_data);

template<typename dreal>
void jamesonrk(
	const struct Flow_options &flo_opts,
    const std::vector<dreal> &area,
    const std::vector<double> &dx,
    class Flow_data<dreal>* const flow_data);

template<typename dreal>
void stepInTime(
	const struct Flow_options &flo_opts,
    const std::vector<dreal> &area,
    const std::vector<double> &dx,
    class Flow_data<dreal>* const flow_data)
{
    if (flo_opts.time_scheme == 0) {
        EulerExplicitStep(flo_opts, area, dx, flow_data);
    } else if (flo_opts.time_scheme == 1) {
		lusgs(flo_opts, area, dx, flow_data);
    } else if (flo_opts.time_scheme == 2) {
        jamesonrk(flo_opts, area, dx, flow_data);
    } else if (flo_opts.time_scheme == 3) {
        abort();//eulerImplicit(area, dx, dt, flow_data);
    } else if (flo_opts.time_scheme == 4) {
        abort();//crankNicolson(area, dx, dt, flow_data);
    }
}
template void stepInTime(
	const struct Flow_options &flo_opts,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    class Flow_data<double>* const flow_data);
template void stepInTime(
	const struct Flow_options &flo_opts,
    const std::vector<adouble> &area,
    const std::vector<double> &dx,
    class Flow_data<adouble>* const flow_data);

// Euler Explicit
template<typename dreal>
void EulerExplicitStep(
	const struct Flow_options &flo_opts,
    const std::vector<dreal> &area,
    const std::vector<double> &dx,
    class Flow_data<dreal>* const flow_data)
{
    getDomainResi(flo_opts, area, flow_data);

	flow_data->old_residual_norm = flow_data->current_residual_norm;
	flow_data->current_residual_norm = norm2(flow_data->residual);
	evaluate_dt(flo_opts, dx, flow_data);

	int n_elem = flo_opts.n_elem;
    for (int k = 0; k < 3; k++){
		for (int i = 1; i < n_elem+1; i++) {
			const int ki_res = i * 3 + k;
			flow_data->W[ki_res] = flow_data->W[ki_res] - (flow_data->dt[i] / dx[i]) * flow_data->residual[ki_res];
		}
	}
    return;
}


// Jameson's 4th order Runge - Kutta Stepping Scheme
template<typename dreal>
void jamesonrk(
	const struct Flow_options &flo_opts,
    const std::vector<dreal> &area,
    const std::vector<double> &dx,
    class Flow_data<dreal>* const flow_data)
{
	int n_elem = flo_opts.n_elem;
    // Initialize First Stage
    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < n_elem+2; i++) {
            const int ki = i * 3 + k;
            flow_data->W_stage[ki] = flow_data->W[ki];
        }
    }
    // 1-4 Stage
    for (int r = 1; r < 5; r++) {
        // Calculate Residuals
		getDomainResi(flo_opts, area, flow_data);
		if (r == 1) {
			flow_data->old_residual_norm = flow_data->current_residual_norm;
			flow_data->current_residual_norm = norm2(flow_data->residual);
			evaluate_dt(flo_opts, dx, flow_data);
		}

        // Step in RK time
        for (int k = 0; k < 3; k++) {
            for (int i = 1; i < n_elem+1; i++) {
                const int ki = i * 3 + k;
                flow_data->W_stage[ki] = flow_data->W[ki] - (flow_data->dt[i-1] / (5 - r)) * flow_data->residual[ki] / dx[i-1];
            }
        }
    }

    for (int k = 0; k < 3; k++) {
        for (int i = 1; i < n_elem+1; i++) {
            const int ki = i * 3 + k;
            flow_data->residual[ki] = (flow_data->W_stage[ki] - flow_data->W[ki]) * dx[i] / flow_data->dt[i];
            flow_data->W[ki] = flow_data->W_stage[ki];
        }
    }
}

template<typename dreal>
void lusgs(
	const struct Flow_options &flo_opts,
    const std::vector<dreal> &area,
    const std::vector<double> &dx,
    class Flow_data<dreal>* const flow_data)
{
	int n_elem = flo_opts.n_elem;

	const double gam = flo_opts.gam;

	using Matrix3 = Eigen::Matrix<dreal, 3, 3>;
	using Vector3 = Eigen::Matrix<dreal, 3, 1>;

	getDomainResi(flo_opts, area, flow_data);

	flow_data->old_residual_norm = flow_data->current_residual_norm;
	flow_data->current_residual_norm = norm2(flow_data->residual);
	evaluate_dt(flo_opts, dx, flow_data);

	Vector3 rhs;
	for (int sweeps = 0; sweeps < 2; sweeps++) {
		for (int row = 1; row < n_elem+1; row++) {

			const dreal u_i = flow_data->W[row*3+1] / flow_data->W[row*3+0];
			const dreal c_i = get_c(gam,flow_data->W[row*3+0],flow_data->W[row*3+1],flow_data->W[row*3+2]);

			const int i_w_p = row+1;
			const int i_w_n = row-1;
			const dreal u_p = flow_data->W[i_w_p*3+1] / flow_data->W[i_w_p*3+0];
			const dreal c_p = get_c(gam,flow_data->W[i_w_p*3+0],flow_data->W[i_w_p*3+1],flow_data->W[i_w_p*3+2]);
			const dreal lambda_p = (u_p+c_p + u_i+c_i)/2.0;

			const dreal u_n = flow_data->W[i_w_n*3+1] / flow_data->W[i_w_n*3+0];
			const dreal c_n = get_c(gam,flow_data->W[i_w_n*3+0],flow_data->W[i_w_n*3+1],flow_data->W[i_w_n*3+2]);
			const dreal lambda_n = (u_n+c_n + u_i+c_i)/2.0;


			const int i_face_p = row;
			const int i_face_n = row-1;
			const dreal area_p = area[i_face_p];
			const dreal area_n = area[i_face_n];

			Vector3 W1;
			W1(0) = flow_data->W[row*3+0];
			W1(1) = flow_data->W[row*3+1];
			W1(2) = flow_data->W[row*3+2];
			Matrix3 jacobian_diag = 0.5*(area_p-area_n)*analytical_flux_jacobian(gam, W1);

			dreal diag_identity = dx[row]/flow_data->dt[row] + 0.5*(lambda_p*area_p + lambda_n*area_n);
			jacobian_diag = jacobian_diag + diag_identity*Matrix3::Identity();


			// Forward sweep
			for (int i_state = 0; i_state < 3; i_state++) {
				rhs(i_state) = -flow_data->residual[row*3+i_state];
			}
			//for (int col = 1; col < row; col++) {
			//  	if (row-col>1) continue;
			for (int col = row-1; col < row; col++) {
				if (col < 1) continue;

				const dreal normal = -1.0;
				W1(0) = flow_data->W[col*3+0] + flow_data->W_stage2[col*3+0];
				W1(1) = flow_data->W[col*3+1] + flow_data->W_stage2[col*3+1];
				W1(2) = flow_data->W[col*3+2] + flow_data->W_stage2[col*3+2];

				rhs = rhs - 0.5*area[col]*normal*WtoF(gam,W1);

				W1(0) = flow_data->W[col*3+0];
				W1(1) = flow_data->W[col*3+1];
				W1(2) = flow_data->W[col*3+2];

				rhs = rhs + 0.5*area[col]*normal*WtoF(gam,W1);

				const dreal u_j = W1(1)/W1(0);
				const dreal c_j = get_c(gam,W1(0),W1(1),W1(2));
				const dreal lambda = (u_i+c_i + u_j+c_j)/2.0;

				for (int i_state = 0; i_state < 3; i_state++) {
					rhs(i_state) = rhs(i_state) + 0.5*lambda*flow_data->W_stage[col*3+i_state]*area[col];
				}
			}
			Vector3 dW = jacobian_diag.fullPivLu().solve(rhs);
			for (int i_state = 0; i_state < 3; i_state++) {
				flow_data->W_stage[row*3+i_state] = dW[i_state];
			}
		}
		for (int row = n_elem; row > 0; row--) {

			dreal w1 = flow_data->W[row*3+0];
			dreal w2 = flow_data->W[row*3+1];
			dreal w3 = flow_data->W[row*3+2];
			Vector3 W1 = VectorToEigen3(w1,w2,w3);

			const dreal u_i = w2 / w1;
			const dreal c_i = get_c(gam,w1,w2,w3);

			const int i_w_p = row+1;
			const int i_w_n = row-1;
			const dreal u_p = flow_data->W[i_w_p*3+1] / flow_data->W[i_w_p*3+0];
			const dreal c_p = get_c(gam,flow_data->W[i_w_p*3+0],flow_data->W[i_w_p*3+1],flow_data->W[i_w_p*3+2]);
			const dreal lambda_p = (u_p+c_p + u_i+c_i)/2.0;

			const dreal u_n = flow_data->W[i_w_n*3+1] / flow_data->W[i_w_n*3+0];
			const dreal c_n = get_c(gam,flow_data->W[i_w_n*3+0],flow_data->W[i_w_n*3+1],flow_data->W[i_w_n*3+2]);
			const dreal lambda_n = (u_n+c_n + u_i+c_i)/2.0;


			const int i_face_p = row;
			const int i_face_n = row-1;
			const dreal area_p = area[i_face_p];
			const dreal area_n = area[i_face_n];
			Matrix3 jacobian_diag = 0.5*(area_p-area_n)*analytical_flux_jacobian(gam, W1);

			dreal diag_identity = dx[row]/flow_data->dt[row] + 0.5*(lambda_p*area_p + lambda_n*area_n);
			jacobian_diag = jacobian_diag + diag_identity*Matrix3::Identity();

			// Backward sweep
			rhs.setZero();
			for (int col = row+1; col < row+2; col++) {
				if (col > n_elem) continue;
			//for (int col = row+1; col < n_elem+1; col++) {
			//  	if (col-row>1) continue;

				const dreal normal = 1.0;
				W1(0) = flow_data->W[col*3+0] + flow_data->W_stage2[col*3+0];
				W1(1) = flow_data->W[col*3+1] + flow_data->W_stage2[col*3+1];
				W1(2) = flow_data->W[col*3+2] + flow_data->W_stage2[col*3+2];

				rhs = rhs - 0.5*area[col]*normal*WtoF(gam,W1);

				W1(0) = flow_data->W[col*3+0];
				W1(1) = flow_data->W[col*3+1];
				W1(2) = flow_data->W[col*3+2];

				rhs = rhs + 0.5*area[col]*normal*WtoF(gam,W1);

				const dreal u_j = W1(1)/W1(0);
				const dreal c_j = get_c(gam,W1(0),W1(1),W1(2));
				const dreal lambda = (u_i+c_i + u_j+c_j)/2.0;

				for (int i_state = 0; i_state < 3; i_state++) {
					rhs(i_state) = rhs(i_state) + 0.5*lambda*flow_data->W_stage2[col*3+i_state]*area[col];
				}
			}
			Vector3 dW = jacobian_diag.fullPivLu().solve(rhs);
			for (int i_state = 0; i_state < 3; i_state++) {
				flow_data->W_stage2[row*3+i_state] = flow_data->W_stage[row*3+i_state] - dW[i_state];
			}
			for (int i_state = 0; i_state < 3; i_state++) {
				flow_data->W[row*3+i_state] = flow_data->W[row*3+i_state] + flow_data->W_stage2[row*3+i_state];
				//std::cout<<flow_data->W[row*3+i_state]<<std::endl;
			}
		} // Row loop
	}
}
template void lusgs( const struct Flow_options &flo_opts, const std::vector<double> &area, const std::vector<double> &dx, class Flow_data<double>* const flow_data);
template void lusgs( const struct Flow_options &flo_opts, const std::vector<adouble> &area, const std::vector<double> &dx, class Flow_data<adouble>* const flow_data);

//void dFdW(
//    std::vector<dreal> &J,
//    dreal Sp, dreal Sm,
//    dreal rho, dreal u, dreal c)
//{
//    J[0] = 0.0;
//    J[1] = 1.0;
//    J[2] = 0.0;
//    J[3] = u * u * (gam - 3.0) / 2.0;
//    J[4] = u * (3.0 - gam);
//    J[5] = gam - 1.0;
//    J[6] = ( pow(u, 3) * (gam - 1.0) 
//           * (gam - 2.0) - 2.0 * u * c * c ) 
//           / (2.0 * (gam - 1.0));
//    J[7] = ( 2.0 * c * c + u * u 
//           * ( -2.0 * gam * gam + 5.0 * gam - 3.0 ) )
//           / (2.0 * (gam - 1.0));
//    J[8] = u * gam;
//
//    dreal Sa = (Sp + Sm) / 2.0;
//    for (int i = 0; i < 9; i++)
//    {
//        J[i] *= Sa;
//    }
//    
//    dreal dpdw[3];
//    dpdw[0] = (gam - 1) / 2.0 * u * u;
//    dpdw[1] = - (gam - 1) * u;
//    dpdw[2] = (gam - 1);
//
//    J[3] -= dpdw[0] * (Sp - Sm);
//    J[4] -= dpdw[1] * (Sp - Sm);
//    J[5] -= dpdw[2] * (Sp - Sm);
//}

//void eulerImplicit(
//    const std::vector<dreal> &area,
//    const std::vector<dreal> &dx,
//    const std::vector<dreal> &dt,
//    Flow_data &flow_data)
//{
//    // Get flow_data.fluxes
//    getFlux(flow_data.fluxes, flow_data.W);
//    // Calculate Residuals
//    getDomainResi(flow_data, area);
//    Eigen::VectorXd RHS(3 * n_elem);
//    RHS.setZero();
//    Eigen::SparseMatrix<dreal> A(3 * n_elem, 3 * n_elem);
//    A = evaldRdW(flow_data.W, dx, dt, area);
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 3; j++) {
//            A.coeffRef(i,j) = 0.0;
//        }
//    }
//    for (int i = 0; i < 3 * n_elem; i++) {
//        A.coeffRef(i,i) += 1.0;
//    }
//
//    int Wik;
//    for (int Wi = 1; Wi < n_elem - 1; Wi++) {
//        for (int Wk = 0; Wk < 3; Wk++) {
//            Wik = Wi * 3 + Wk;
//            RHS[Wik] += - dt[Wi] / dx[Wi] * flow_data.residual[Wik];
//
//            if (Wi == 0 || Wi == n_elem-1) {
//                RHS[Wik] = -dt[Wi] / dx[Wi] * flow_data.residual[Wik];
//            }
//        }
//    }
//    VectorXd Wt(3 * n_elem);
//    SparseLU <SparseMatrix<dreal>, COLAMDOrdering< int > > slusolver1;
//    slusolver1.compute(A);
//    if (slusolver1.info() != 0)
//        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
//    Wt = slusolver1.solve(RHS);
//    dreal currentR = 0;
//    for (int i = 0; i < n_elem; i++)
//        currentR += flow_data.residual[i * 3 + 0] * flow_data.residual[i * 3 + 0];
//    currentR = sqrt(currentR);
//    dreal alpha = 1;
//    if (1.0/currentR > 1e2)
//    {
//        alpha = 1.25 * pow(log10(1.0/currentR/1e1), 3);
//        std::cout<<alpha<<std::endl;
//    }
//    std::cout<<alpha<<std::endl;
//    for (int k = 0; k < 3; k++) {
//        for (int i = 1; i < n_elem - 1; i++) {
//            ki = i * 3 + k;
//            flow_data.W[ki] += Wt[ki] * alpha;
//            flow_data.residual[ki] = Wt[ki] * dx[i] / dt[i];
//        }
//    }
//}
//void crankNicolson(
//    const std::vector<dreal> &area,
//    const std::vector<dreal> &dx,
//    const std::vector<dreal> &dt,
//    Flow_data &flow_data)
//{
//    // Get flow_data.fluxes
//    getFlux(flow_data.fluxes, flow_data.W);
//    // Calculate Residuals
//    getDomainResi(flow_data, area);
//    Eigen::VectorXd RHS(3 * n_elem);
//    RHS.setZero();
//    Eigen::SparseMatrix<dreal> A(3 * n_elem, 3 * n_elem);
//    A = 0.5 * evaldRdW(flow_data.W, dx, dt, area);
//    for (int i = 0; i < 3; i++) {
//        for (int j = 0; j < 3; j++) {
//            A.coeffRef(i,j) = 0.0;
//        }
//    }
//    for (int i = 0; i < 3 * n_elem; i++) {
//        A.coeffRef(i,i) += 1.0;
//    }
//
//    int Wik;
//    for (int Wi = 1; Wi < n_elem - 1; Wi++) {
//        for (int Wk = 0; Wk < 3; Wk++) {
//            Wik = Wi * 3 + Wk;
//            RHS[Wik] += - dt[Wi] / dx[Wi] * flow_data.residual[Wik];
//
//            if (Wi == 0 || Wi == n_elem-1) {
//                RHS[Wik] = -dt[Wi] / dx[Wi] * flow_data.residual[Wik];
//            }
//        }
//    }
////  std::cout<<A<<std::endl;
////  std::cout<<RHS<<std::endl;
////  *((unsigned int*)0) = 0xDEAD;
//    VectorXd Wt(3 * n_elem);
////  Wt = solveGMRES(A, RHS);
//    SparseLU <SparseMatrix<dreal>, COLAMDOrdering< int > > slusolver1;
//    slusolver1.compute(A);
//    if (slusolver1.info() != 0)
//        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
//    Wt = slusolver1.solve(RHS);
//    dreal currentR = 0;
//    for (int i = 0; i < n_elem; i++)
//        currentR += flow_data.residual[i * 3 + 0] * flow_data.residual[i * 3 + 0];
//    currentR = sqrt(currentR);
//    dreal alpha = 1;
//    if (1.0/currentR > 1e2) {
//        alpha = 1.1 * pow(log10(1.0/currentR/1e1), 3.0);
//        std::cout<<alpha<<std::endl;
//    }
//    std::cout<<alpha<<std::endl;
//    for (int k = 0; k < 3; k++) {
//        for (int i = 1; i < n_elem - 1; i++) {
//            ki = i * 3 + k;
//            flow_data.W[ki] += Wt[ki] * alpha;
//            flow_data.residual[ki] = Wt[ki] * dx[i] / dt[i];
//        }
//    }
//}

//void rk4(
//    const std::vector<dreal> &area,
//    const std::vector<dreal> &dx,
//    const std::vector<dreal> &dt,
//    Flow_data &flow_data);
//
// 4th order Runge - Kutta Stepping Scheme
// Need to figure out how to store more W
//void rk4(
//    const std::vector<dreal> &area,
//    const std::vector<dreal> &dx,
//    const std::vector<dreal> &dt,
//    Flow_data &flow_data)
//{
//    // Residual 0
//    getFlux(flow_data.fluxes, flow_data.W);
//    getDomainResi(flow_data.W, flow_data.fluxes, area, Resi0);
//    // RK1
//    for (int k = 0; k < 3; k++)
//    {
//        for (int i = 1; i < n_elem - 1; i++)
//        {
//            ki = i * 3 + k;
//            W1[ki] = flow_data.W[ki] - (dt[i] / 2) * Resi0[ki] / dx[i];
//        }
//        W1[0 * 3 + k] = flow_data.W[0 * 3 + k];
//        W1[(n_elem - 1) * 3 + k] = flow_data.W[(n_elem - 1) * 3 + k];
//    }
//
//    // Residual 1
//    getFlux(flow_data.fluxes, W1);
//    getDomainResi(W1, flow_data.fluxes, area, Resi1);
//
//    // RK2
//    for (int k = 0; k < 3; k++)
//    {
//        for (int i = 1; i < n_elem - 1; i++)
//        {
//            ki = i * 3 + k;
//            W2[ki] = flow_data.W[ki] - (dt[i] / 2) * Resi1[ki] / dx[i];
//        }
//        W2[0 * 3 + k] = flow_data.W[0 * 3 + k];
//        W2[(n_elem - 1) * 3 + k] = flow_data.W[(n_elem - 1) * 3 + k];
//    }
//
//    // Residual 2
//    getFlux(flow_data.fluxes, W2);
//    getDomainResi(W2, flow_data.fluxes, area, Resi2);
//
//    // RK3
//    for (int k = 0; k < 3; k++)
//    {
//        for (int i = 1; i < n_elem - 1; i++)
//        {
//            ki = i * 3 + k;
//            W3[ki] = flow_data.W[ki] - (dt[i] / 2) * Resi2[ki] / dx[i];
//        }
//    }
//
//    for (int k = 0; k < 3; k++)
//    {
//        for (int i = 1; i < n_elem - 1; i++)
//        {
//            ki = i * 3 + k;
//            Wtemp[ki] = ((dreal)1.0 / 6.0) * (flow_data.W[ki] + 2 * W1[ki] + 2 * W2[ki] + W3[ki]);
//            //flow_data.residual[ki] = (2 * Resi0[ki] + 2 * Resi1[ki] + Resi2[ki]) / 6.0;
//            flow_data.residual[ki] = (Wtemp[ki] - flow_data.W[ki]) * dx[i] / dt[i];
//            flow_data.W[ki] = Wtemp[ki];
//        }
//    }
//}

