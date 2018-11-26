// Time Stepping Schemes
#include "timestep.hpp"
#include "structures.hpp"
#include<vector>
#include<math.h>
#include<iostream>
#include<Eigen/Core>
#include<Eigen/LU>
#include<Eigen/Sparse>
#include<Eigen/Dense> // eigenvalues
#include "flux.hpp"
#include "convert.hpp"
//#include "residuald1.hpp"
#include "boundary_conditions.hpp"
#include "boundary_gradient2.hpp"
#include "derivatives.hpp"
#include<adolc/adolc.h>
#include"adolc_eigen.hpp"
#include"residual_derivatives.hpp"
#include"residuald1.hpp"
#include"solve_linear.hpp"

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
void BackwardEuler(
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
    //if (flow_data->current_residual_norm > 1e-5) {
    //    BackwardEuler(flo_opts, area, dx, flow_data);//eulerImplicit(area, dx, dt, flow_data);
    //} else {
	//	lusgs(flo_opts, area, dx, flow_data);
    //}
    //return;
    if (flo_opts.time_scheme == 0) {
        EulerExplicitStep(flo_opts, area, dx, flow_data);
    } else if (flo_opts.time_scheme == 1) {
		lusgs(flo_opts, area, dx, flow_data);
    } else if (flo_opts.time_scheme == 2) {
        jamesonrk(flo_opts, area, dx, flow_data);
    } else if (flo_opts.time_scheme == 3) {
        BackwardEuler(flo_opts, area, dx, flow_data);//eulerImplicit(area, dx, dt, flow_data);
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

//void stepInTime2(
//	const struct Flow_options &flo_opts,
//    const std::vector<double> &area,
//    const std::vector<double> &dx,
//    class Flow_data<double>* const flow_data)
//{
//	int n_elem = flo_opts.n_elem;
//
//	const double gam = flo_opts.gam;
//	const double half = 0.5;
//
//	using Matrix3<dreal> = Eigen::Matrix3<double, 3, 3>;
//	using Vector3<dreal> = Eigen::Matrix3<double, 3, 1>;
//	Matrix3<dreal> identity = Matrix3<dreal>::Identity();
//
//	Eigen::SparseMatrix<double> dRdW = eval_dRdW_dRdX_adolc(flo_opts, area_double, flow_data_double);
//	if (print_jac) {
//				std::cout<<dRdW;
//				std::cout<<dRdW<<std::endl;
//	}
//
//	getDomainResi(flo_opts, area, flow_data);
//
//	flow_data->old_residual_norm = flow_data->current_residual_norm;
//	flow_data->current_residual_norm = norm2(flow_data->residual);
//	evaluate_dt(flo_opts, dx, flow_data);
//}

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

template<>
void BackwardEuler(
	const struct Flow_options &flo_opts,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    class Flow_data<double>* const flow_data)
{
	int n_elem = flo_opts.n_elem;

    getDomainResi(flo_opts, area, flow_data);
	flow_data->old_residual_norm = flow_data->current_residual_norm;
	flow_data->current_residual_norm = norm2(flow_data->residual);
	evaluate_dt(flo_opts, dx, flow_data);

	//Eigen::SparseMatrix<double> dRdW = eval_dRdW_dRdX_adolc(flo_opts, area, *flow_data);
	Eigen::SparseMatrix<double> dRdW = evaldRdW_FD(area, flo_opts, *flow_data);
    Eigen::VectorXd minusR(n_elem*3);
    for (int k = 0; k < 3; k++){
		for (int i = 1; i < n_elem+1; i++) {
			const int ki = i * 3 + k;
			const int kim = (i-1) * 3 + k;
			minusR(kim) = -flow_data->residual[ki];
            dRdW.coeffRef(kim,kim) += dx[i]/flow_data->dt[i];
		}
	}
    Eigen::VectorXd dW = solve_linear(dRdW, minusR, 0, 1e-9);

    for (int k = 0; k < 3; k++){
		for (int i = 1; i < n_elem+1; i++) {
			const int ki = i * 3 + k;
			const int kim = (i-1) * 3 + k;
			flow_data->W[ki] = flow_data->W[ki] + dW[kim];
		}
	}
    return;
}
void BackwardEuler(
	const struct Flow_options &flo_opts,
    const std::vector<adouble> &area,
    const std::vector<double> &dx,
    class Flow_data<adouble>* const flow_data)
{return;}


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
            flow_data->W_stage1[ki] = flow_data->W[ki];
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
                flow_data->W_stage1[ki] = flow_data->W[ki] - (flow_data->dt[i-1] / (5 - r)) * flow_data->residual[ki] / dx[i-1];
            }
        }
    }

    for (int k = 0; k < 3; k++) {
        for (int i = 1; i < n_elem+1; i++) {
            const int ki = i * 3 + k;
            flow_data->residual[ki] = (flow_data->W_stage1[ki] - flow_data->W[ki]) * dx[i] / flow_data->dt[i];
            flow_data->W[ki] = flow_data->W_stage1[ki];
        }
    }
}

double valueof(double a) { return a; }
double valueof(adouble a) { return a.value(); }

template<typename dreal>
void lusgs(
	const struct Flow_options &flo_opts,
    const std::vector<dreal> &area,
    const std::vector<double> &dx,
    class Flow_data<dreal>* const flow_data)
{
	int n_elem = flo_opts.n_elem;

	const double gam = flo_opts.gam;
	const dreal half = 0.5;

	Matrix3<dreal> identity = Matrix3<dreal>::Identity();

	bool print_jac = false;// print_jac=true;
	std::vector<double> area_double(n_elem+1);
	Flow_data<double> flow_data_double(n_elem);
	for (int row = 0; row < n_elem+2; row++) {
		for (int i_state = 0; i_state < 3; i_state++) {
			flow_data_double.W[row*3+i_state] = valueof(flow_data->W[row*3+i_state]);
		}
	}
	for (int row = 0; row < n_elem+1; row++) {
		area_double[row] = valueof(area[row]);
	}
	Eigen::SparseMatrix<double> dRdW = eval_dRdW_dRdX_adolc(flo_opts, area_double, flow_data_double);
	
	getDomainResi(flo_opts, area, flow_data);

	flow_data->old_residual_norm = flow_data->current_residual_norm;
	flow_data->current_residual_norm = norm2(flow_data->residual);
	evaluate_dt(flo_opts, dx, flow_data);

    VectorXd dW_total(3*n_elem);
    VectorXd residual(3*n_elem);

	Vector3<dreal> rhs;
	for (int row = 1; row < n_elem+1; row++) {
		for (int i_state = 0; i_state < 3; i_state++) {
			flow_data->W_stage1[row*3+i_state] = 0.0; // dW_forward
			flow_data->W_stage2[row*3+i_state] = 0.0; // dW_backward

            residual((row-1)*3+i_state) = flow_data->residual[row*3+i_state];
		}
	}
	const int n_sweeps = 20;
	int print_row = -1;
	for (int i_sweep = 0; i_sweep < n_sweeps; i_sweep++) {

		dreal normdw = 0;
		// Forward sweep
		for (int row = 1; row < n_elem+1; row++) {

			const int i_w_p = row+1;
			const int i_w_n = row-1;

			const int i_face_p = row;
			const int i_face_n = row-1;
			const dreal area_p = area[i_face_p];
			const dreal area_n = area[i_face_n];

			Vector3<dreal> Wd, Wn, Wp, dW_stage;
			Wd(0) = flow_data->W[row*3+0];
			Wd(1) = flow_data->W[row*3+1];
			Wd(2) = flow_data->W[row*3+2];
			Wn(0) = flow_data->W[i_w_n*3+0];
			Wn(1) = flow_data->W[i_w_n*3+1];
			Wn(2) = flow_data->W[i_w_n*3+2];
			Wp(0) = flow_data->W[i_w_p*3+0];
			Wp(1) = flow_data->W[i_w_p*3+1];
			Wp(2) = flow_data->W[i_w_p*3+2];
			Matrix3<dreal> jacobian_diag = (analytical_flux_jacobian(gam, Wd))*half;
			jacobian_diag(1,0) = jacobian_diag(1,0) - get_dpdw1(gam, Wd(0), Wd(1), Wd(2));
			jacobian_diag(1,1) = jacobian_diag(1,1) - get_dpdw2(gam, Wd(0), Wd(1), Wd(2));
			jacobian_diag(1,2) = jacobian_diag(1,2) - get_dpdw3(gam, Wd(0), Wd(1), Wd(2));
			jacobian_diag = (area_p-area_n)*jacobian_diag;
			
			Matrix3<dreal> d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, Wn(0), Wn(1), Wn(2), Wd(0), Wd(1), Wd(2));
			jacobian_diag -= half * flo_opts.scalar_d_eps * (-d_lambdaw_dw) * area_n;

			d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, Wd(0), Wd(1), Wd(2), Wp(0), Wp(1), Wp(2));
			jacobian_diag += half * flo_opts.scalar_d_eps * (-d_lambdaw_dw) * area_p;

			if (row == 1) { // Add boundary contribution
				Matrix3<dreal> dWndW = inletBC_gradient(flo_opts, flow_data);
				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, Wn(0), Wn(1), Wn(2), Wd(0), Wd(1), Wd(2));
				d_lambdaw_dw *= flo_opts.scalar_d_eps;
				Matrix3<dreal> dRdWn = -area_n * half * (analytical_flux_jacobian(gam, Wn) - d_lambdaw_dw);
				if(row==print_row) std::cout<<"dWndW "<<std::endl<<dWndW<<std::endl;
				if(row==print_row) std::cout<<"dRdWn "<<std::endl<<dRdWn<<std::endl;
				jacobian_diag = jacobian_diag + dRdWn*dWndW;
			}
			if (row == n_elem) { // Add boundary contribution
				Matrix3<dreal> dWpdW = outletBC_gradient(flo_opts, flow_data);
				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, Wd(0), Wd(1), Wd(2), Wp(0), Wp(1), Wp(2));
				d_lambdaw_dw *= flo_opts.scalar_d_eps;
				Matrix3<dreal> dRdWp = area_p * half * (analytical_flux_jacobian(gam, Wp) - d_lambdaw_dw);
				if(row==print_row) std::cout<<"dWpdW "<<std::endl<<dWpdW<<std::endl;
				if(row==print_row) std::cout<<"dRdWp "<<std::endl<<dRdWp<<std::endl;
				jacobian_diag = jacobian_diag + dRdWp*dWpdW;
			}
			Matrix3<dreal> Middle;
			for (int i = 0; i<3; i++) {
				for (int j = 0; j<3; j++) {
					Middle(i,j) = dRdW.coeff((row-1)*3+i,(row-1)*3+j);
				}
			}
			if(row==print_row) std::cout<<"Row "<<row<<std::endl<<jacobian_diag<<std::endl;

			if(print_jac) std::cout<<"Middle"<<std::endl;
			if(print_jac) std::cout<<row<<std::endl<<jacobian_diag<<std::endl;
			if(print_jac) std::cout<<Middle<<std::endl<<std::endl;
			if(print_jac) std::cout<<jacobian_diag-Middle<<std::endl<<std::endl;
			if(print_jac) if((jacobian_diag-Middle).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;


			// Add pseudo-timestep
			dreal diag_identity = dx[row]/flow_data->dt[row];
			jacobian_diag = jacobian_diag + diag_identity*identity;
			if(row==print_row) std::cout<<"Row "<<row<<std::endl<<jacobian_diag<<std::endl;


			for (int i_state = 0; i_state < 3; i_state++) {
				rhs(i_state) = -flow_data->residual[row*3+i_state];
			}
			if(row==print_row) std::cout<<"rhs "<<rhs<<std::endl<<std::endl;
			for (int lower_index = row-1; lower_index < row; lower_index++) {
				// U*dW
				if (row!=n_elem) {
					int i_p = row+1;
					d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, Wd(0), Wd(1), Wd(2), Wp(0), Wp(1), Wp(2));
					d_lambdaw_dw *= flo_opts.scalar_d_eps;
					Matrix3<dreal> dRdWp = area_p * half * (analytical_flux_jacobian(gam, Wp) - d_lambdaw_dw);

					dW_stage(0) = flow_data->W_stage2[i_p*3+0];
					dW_stage(1) = flow_data->W_stage2[i_p*3+1];
					dW_stage(2) = flow_data->W_stage2[i_p*3+2];

					Matrix3<dreal> Upper;
					for (int i = 0; i<3; i++) {
						for (int j = 0; j<3; j++) {
							Upper(i,j) = dRdW.coeff((row-1)*3+i,(i_p-1)*3+j);
						}
					}
					if(print_jac) std::cout<<"Upper"<<std::endl;
					if(print_jac) std::cout<<dRdWp<<std::endl;
					if(print_jac) std::cout<<Upper<<std::endl<<std::endl;
					if(print_jac) std::cout<<dRdWp-Upper<<std::endl<<std::endl;
					if(print_jac) if((dRdWp-Upper).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;


					rhs -= dRdWp*dW_stage;
				}

				if (lower_index < 1) continue;

				// -L*dW
				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, Wn(0), Wn(1), Wn(2), Wd(0), Wd(1), Wd(2));
				d_lambdaw_dw *= flo_opts.scalar_d_eps;
				Matrix3<dreal> dRdWn = -area_n * half * (analytical_flux_jacobian(gam, Wn) - d_lambdaw_dw);

				dW_stage(0) = flow_data->W_stage1[lower_index*3+0];
				dW_stage(1) = flow_data->W_stage1[lower_index*3+1];
				dW_stage(2) = flow_data->W_stage1[lower_index*3+2];

				rhs -= dRdWn*dW_stage;

				
				Matrix3<dreal> Lower;
				for (int i = 0; i<3; i++) {
					for (int j = 0; j<3; j++) {
						Lower(i,j) = dRdW.coeff((row-1)*3+i,(lower_index-1)*3+j);
					}
				}
				if(print_jac) std::cout<<"Lower"<<std::endl;
				if(print_jac) std::cout<<dRdWn<<std::endl;
				if(print_jac) std::cout<<Lower<<std::endl<<std::endl;
				if(print_jac) std::cout<<dRdWn-Lower<<std::endl<<std::endl;
				if(print_jac) if((dRdWn-Lower).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;

			}
			Vector3<dreal> dW = jacobian_diag.fullPivLu().solve(rhs);
			for (int i_state = 0; i_state < 3; i_state++) {
				flow_data->W_stage1[row*3+i_state] = dW[i_state];
				normdw += dW.dot(dW);
                dW_total((row-1)*3+i_state) = dW[i_state];
			}
			if(row==print_row) std::cout<<"dW* "<<dW<<std::endl<<std::endl;
			if(row==print_row) std::cout<<"A "<<jacobian_diag<<std::endl<<std::endl;
			if(row==print_row) std::cout<<"rhs "<<rhs<<std::endl<<std::endl;
			//if(row==print_row) abort();
		} // Row loop
		std::cout<<"forward i_sweep"<<i_sweep<<" normdw "<<normdw<<std::endl;
		std::cout<<"backward residual"<<(dRdW*dW_total+residual).norm()<<std::endl;

		// Backward sweep
		normdw = 0;
		for (int row = n_elem; row > 0; row--) {

			const int i_w_p = row+1;
			const int i_w_n = row-1;

			const int i_face_p = row;
			const int i_face_n = row-1;
			const dreal area_p = area[i_face_p];
			const dreal area_n = area[i_face_n];

			Vector3<dreal> Wd, Wn, Wp, dW_stage;
			Wd(0) = flow_data->W[row*3+0];
			Wd(1) = flow_data->W[row*3+1];
			Wd(2) = flow_data->W[row*3+2];
			Wn(0) = flow_data->W[i_w_n*3+0];
			Wn(1) = flow_data->W[i_w_n*3+1];
			Wn(2) = flow_data->W[i_w_n*3+2];
			Wp(0) = flow_data->W[i_w_p*3+0];
			Wp(1) = flow_data->W[i_w_p*3+1];
			Wp(2) = flow_data->W[i_w_p*3+2];
			Matrix3<dreal> jacobian_diag = (analytical_flux_jacobian(gam, Wd))*half;
			jacobian_diag(1,0) = jacobian_diag(1,0) - get_dpdw1(gam, Wd(0), Wd(1), Wd(2));
			jacobian_diag(1,1) = jacobian_diag(1,1) - get_dpdw2(gam, Wd(0), Wd(1), Wd(2));
			jacobian_diag(1,2) = jacobian_diag(1,2) - get_dpdw3(gam, Wd(0), Wd(1), Wd(2));
			jacobian_diag = (area_p-area_n)*jacobian_diag;
			
			Matrix3<dreal> d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, Wn(0), Wn(1), Wn(2), Wd(0), Wd(1), Wd(2));
			jacobian_diag -= half * flo_opts.scalar_d_eps * (-d_lambdaw_dw) * area_n;

			d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, Wd(0), Wd(1), Wd(2), Wp(0), Wp(1), Wp(2));
			jacobian_diag += half * flo_opts.scalar_d_eps * (-d_lambdaw_dw) * area_p;

			if (row == 1) { // Add boundary contribution
				Matrix3<dreal> dWndW = inletBC_gradient(flo_opts, flow_data);
				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, Wn(0), Wn(1), Wn(2), Wd(0), Wd(1), Wd(2));
				d_lambdaw_dw *= flo_opts.scalar_d_eps;
				Matrix3<dreal> dRdWn = -area_n * half * (analytical_flux_jacobian(gam, Wn) - d_lambdaw_dw);
				if(row==print_row) std::cout<<"dWndW "<<std::endl<<dWndW<<std::endl;
				if(row==print_row) std::cout<<"dRdWn "<<std::endl<<dRdWn<<std::endl;
				jacobian_diag = jacobian_diag + dRdWn*dWndW;
			}
			if (row == n_elem) { // Add boundary contribution
				Matrix3<dreal> dWpdW = outletBC_gradient(flo_opts, flow_data);
				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, Wd(0), Wd(1), Wd(2), Wp(0), Wp(1), Wp(2));
				d_lambdaw_dw *= flo_opts.scalar_d_eps;
				Matrix3<dreal> dRdWp = area_p * half * (analytical_flux_jacobian(gam, Wp) - d_lambdaw_dw);
				if(row==print_row) std::cout<<"dWpdW "<<std::endl<<dWpdW<<std::endl;
				if(row==print_row) std::cout<<"dRdWp "<<std::endl<<dRdWp<<std::endl;
				jacobian_diag = jacobian_diag + dRdWp*dWpdW;
			}
			Matrix3<dreal> Middle;
			for (int i = 0; i<3; i++) {
				for (int j = 0; j<3; j++) {
					Middle(i,j) = dRdW.coeff((row-1)*3+i,(row-1)*3+j);
				}
			}
			if(row==print_row) std::cout<<"Row "<<row<<std::endl<<jacobian_diag<<std::endl;

			if(print_jac) std::cout<<"Middle"<<std::endl;
			if(print_jac) std::cout<<row<<std::endl<<jacobian_diag<<std::endl;
			if(print_jac) std::cout<<Middle<<std::endl<<std::endl;
			if(print_jac) std::cout<<jacobian_diag-Middle<<std::endl<<std::endl;
			if(print_jac) if((jacobian_diag-Middle).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;


			// Add pseudo-timestep
			dreal diag_identity = dx[row]/flow_data->dt[row];
			jacobian_diag = jacobian_diag + diag_identity*identity;
			if(row==print_row) std::cout<<"Row "<<row<<std::endl<<jacobian_diag<<std::endl;

			//jacobian_diag = 10*diag_identity*identity;


			for (int i_state = 0; i_state < 3; i_state++) {
				rhs(i_state) = -flow_data->residual[row*3+i_state];
			}
			// Backward sweep
			for (int upper_index = row+1; upper_index < row+2; upper_index++) {
				// -L*dW
				if (row!=1) {
					int i_n = row-1;
					d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, Wn(0), Wn(1), Wn(2), Wd(0), Wd(1), Wd(2));
					d_lambdaw_dw *= flo_opts.scalar_d_eps;
					Matrix3<dreal> dRdWn = -area_n * half * (analytical_flux_jacobian(gam, Wn) - d_lambdaw_dw);

					Matrix3<dreal> Lower;
					for (int i = 0; i<3; i++) {
						for (int j = 0; j<3; j++) {
							Lower(i,j) = dRdW.coeff((row-1)*3+i,(i_n-1)*3+j);
						}
					}
					if(print_jac) std::cout<<"Lower"<<std::endl;
					if(print_jac) std::cout<<dRdWn<<std::endl;
					if(print_jac) std::cout<<Lower<<std::endl<<std::endl;
					if(print_jac) std::cout<<dRdWn-Lower<<std::endl<<std::endl;
					if(print_jac) if((dRdWn-Lower).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;

					dW_stage(0) = flow_data->W_stage1[i_n*3+0];
					dW_stage(1) = flow_data->W_stage1[i_n*3+1];
					dW_stage(2) = flow_data->W_stage1[i_n*3+2];

					rhs -= dRdWn*dW_stage;
				}
				if (upper_index > n_elem) continue;

				// -U*dW
				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, Wd(0), Wd(1), Wd(2), Wp(0), Wp(1), Wp(2));
				d_lambdaw_dw *= flo_opts.scalar_d_eps;
				Matrix3<dreal> dRdWp = area_p * half * (analytical_flux_jacobian(gam, Wp) - d_lambdaw_dw);

				Matrix3<dreal> Upper;
				for (int i = 0; i<3; i++) {
					for (int j = 0; j<3; j++) {
						Upper(i,j) = dRdW.coeff((row-1)*3+i,(upper_index-1)*3+j);
					}
				}
				if(print_jac) std::cout<<"Upper"<<std::endl;
				if(print_jac) std::cout<<dRdWp<<std::endl;
				if(print_jac) std::cout<<Upper<<std::endl<<std::endl;
				if(print_jac) std::cout<<dRdWp-Upper<<std::endl<<std::endl;
				if(print_jac) if((dRdWp-Upper).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;

				dW_stage(0) = flow_data->W_stage2[upper_index*3+0];
				dW_stage(1) = flow_data->W_stage2[upper_index*3+1];
				dW_stage(2) = flow_data->W_stage2[upper_index*3+2];

				rhs -= dRdWp*dW_stage;

			}
			Vector3<dreal> dW = jacobian_diag.fullPivLu().solve(rhs);
			for (int i_state = 0; i_state < 3; i_state++) {
				flow_data->W_stage2[row*3+i_state] = dW[i_state];
				normdw += dW.dot(dW);
                dW_total((row-1)*3+i_state) = dW[i_state];
			}
		} // Row loop
		std::cout<<"backward i_sweep"<<i_sweep<<" normdw "<<sqrt(normdw)<<std::endl;
		std::cout<<"backward residual"<<(dRdW*dW_total+residual).norm()<<std::endl;
		if(print_jac) abort();
	} // Sweep loop
	for (int row = 1; row < n_elem+1; row++) {
		for (int i_state = 0; i_state < 3; i_state++) {
			flow_data->W[row*3+i_state] = flow_data->W[row*3+i_state] + flow_data->W_stage2[row*3+i_state];
		}
	}
}
template void lusgs( const struct Flow_options &flo_opts, const std::vector<double> &area, const std::vector<double> &dx, class Flow_data<double>* const flow_data);
//template<typename dreal>
//void lusgs(
//	const struct Flow_options &flo_opts,
//    const std::vector<dreal> &area,
//    const std::vector<double> &dx,
//    class Flow_data<dreal>* const flow_data)
//{
//	int n_elem = flo_opts.n_elem;
//
//	const double gam = flo_opts.gam;
//	const dreal half = 0.5;
//
//	using Matrix3<dreal> = Eigen::Matrix3<dreal, 3, 3>;
//	using Vector3<dreal> = Eigen::Matrix3<dreal, 3, 1>;
//	Matrix3<dreal> identity = Matrix3<dreal>::Identity();
//
//	bool print_jac = false;// print_jac=true;
//	std::vector<double> area_double(n_elem+1);
//	Flow_data<double> flow_data_double(n_elem);
//	for (int row = 0; row < n_elem+2; row++) {
//		for (int i_state = 0; i_state < 3; i_state++) {
//			flow_data_double.W[row*3+i_state] = valueof(flow_data->W[row*3+i_state]);
//		}
//	}
//	for (int row = 0; row < n_elem+1; row++) {
//		area_double[row] = valueof(area[row]);
//	}
//	Eigen::SparseMatrix<double> dRdW = eval_dRdW_dRdX_adolc(flo_opts, area_double, flow_data_double);
//	Eigen::MatrixXd dRdW_dense = Eigen::MatrixXd(dRdW);
//	Eigen::MatrixXd L(3*n_elem,3*n_elem), U(3*n_elem,3*n_elem), D(3*n_elem,3*n_elem);
//	L.setZero();
//	U.setZero();
//	D.setZero();
//	
//
//	//std::cout<<dRdW<<std::endl;
//	//for (int row = 0; row < 3*n_elem; row++){
//	//	double diag = valueof(dx[row/3]/flow_data->dt[row/3]);
//	//	diag += abs(dRdW.coeff(row,row));
//	//	for (int col=0; col<3*n_elem;col++) {
//	//		if(row==col) continue;
//	//		diag -= abs(dRdW.coeff(row,col));
//	//	}
//	//	std::cout<<"diag row="<<row<<" dominance "<<diag<<std::endl;
//	//}
//	if (print_jac) {
//				std::cout<<dRdW;
//				std::cout<<dRdW<<std::endl;
//	}
//
//	getDomainResi(flo_opts, area, flow_data);
//
//	flow_data->old_residual_norm = flow_data->current_residual_norm;
//	flow_data->current_residual_norm = norm2(flow_data->residual);
//	evaluate_dt(flo_opts, dx, flow_data);
//
//	for (int row = 0; row < n_elem; row++){
//		double d = valueof(dx[row]/flow_data->dt[row]);
//		for(int i = 0; i<3; i++) {
//			for(int j = 0; j<3; j++) {
//				D(row*3+i,row*3+j) = dRdW_dense(row*3+i,row*3+j) + d;
//				if(row!=0) {
//					int rowm = row-1;
//					U(rowm*3+i,row*3+j) = dRdW_dense(rowm*3+i,row*3+j);
//				}
//				if(row!=n_elem-1) {
//					int rowp = row+1;
//					L(rowp*3+i,row*3+j) = dRdW_dense(rowp*3+i,row*3+j);
//				}
//			}
//		}
//	}
//	//std::cout<<((D+L+U)-dRdW_dense).norm();
//	//Eigen::MatrixXd T_forward= (D+L).inverse()*U;
//	//Eigen::MatrixXd T_backward= (D+U).inverse()*L;
//
//	//Eigen::MatrixXd T = T_backward*T_forward;
//	//std::cout<<"Forward eigenvalues"<<std::endl;
//	//std::cout<<T_forward.eigenvalues()<<std::endl;
//	//std::cout<<"Backward eigenvalues"<<std::endl;
//	//std::cout<<T_backward.eigenvalues()<<std::endl;
//	//std::cout<<"Total eigenvalues"<<std::endl;
//	//std::cout<<T.eigenvalues()<<std::endl;
//
//	Vector3<dreal> rhs;
//	for (int row = 0; row < n_elem+1; row++) {
//		for (int i_state = 0; i_state < 3; i_state++) {
//			flow_data->W_stage1[row*3+i_state] = 0.0;
//			flow_data->W_stage2[row*3+i_state] = 0.0;
//		}
//	}
//	const int n_sweeps = 1;
//	int print_row = -1;
//	for (int i_sweep = 0; i_sweep < n_sweeps; i_sweep++) {
//
//		// Forward sweep
//		for (int row = 1; row < n_elem+1; row++) {
//
//			const dreal u_i = flow_data->W[row*3+1] / flow_data->W[row*3+0];
//			const dreal c_i = get_c(gam,flow_data->W[row*3+0],flow_data->W[row*3+1],flow_data->W[row*3+2]);
//
//			const int i_w_p = row+1;
//			const int i_w_n = row-1;
//			const dreal u_n = flow_data->W[i_w_n*3+1] / flow_data->W[i_w_n*3+0];
//			const dreal c_n = get_c(gam,flow_data->W[i_w_n*3+0],flow_data->W[i_w_n*3+1],flow_data->W[i_w_n*3+2]);
//
//			const dreal lambda_n = (u_n+c_n + u_i+c_i)/2.0;
//
//			const int i_face_p = row;
//			const int i_face_n = row-1;
//			const dreal area_p = area[i_face_p];
//			const dreal area_n = area[i_face_n];
//
//			Vector3<dreal> W1;
//			W1(0) = flow_data->W[row*3+0];
//			W1(1) = flow_data->W[row*3+1];
//			W1(2) = flow_data->W[row*3+2];
//			Matrix3<dreal> jacobian_diag = (analytical_flux_jacobian(gam, W1))*half;
//			jacobian_diag(1,0) = jacobian_diag(1,0) - get_dpdw1(gam, W1(0), W1(1), W1(2));
//			jacobian_diag(1,1) = jacobian_diag(1,1) - get_dpdw2(gam, W1(0), W1(1), W1(2));
//			jacobian_diag(1,2) = jacobian_diag(1,2) - get_dpdw3(gam, W1(0), W1(1), W1(2));
//			jacobian_diag = (area_p-area_n)*jacobian_diag;
//			
//			Matrix3<dreal> d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, flow_data->W[i_w_n*3+0], flow_data->W[i_w_n*3+1], flow_data->W[i_w_n*3+2], W1(0), W1(1), W1(2));
//			jacobian_diag += half * flo_opts.scalar_d_eps * d_lambdaw_dw * area_n;
//
//			d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, W1(0), W1(1), W1(2), flow_data->W[i_w_p*3+0], flow_data->W[i_w_p*3+1], flow_data->W[i_w_p*3+2]);
//			jacobian_diag -= half * flo_opts.scalar_d_eps * d_lambdaw_dw * area_p;
//
//
//
//			if (row == 1) { // Add boundary contribution
//				//Vector3<dreal> W2;
//				//W2(0) = flow_data->W[(row-1)*3+0];
//				//W2(1) = flow_data->W[(row-1)*3+1];
//				//W2(2) = flow_data->W[(row-1)*3+2];
//				//Matrix3<dreal> dWidWd = inletBC_gradient(flo_opts, flow_data);
//				//Matrix3<dreal> dRddWi = -half*(area_n) * (analytical_flux_jacobian(gam, W2) - lambda_n*identity);
//				//jacobian_diag = jacobian_diag + dRddWi*dWidWd;
//
//				Vector3<dreal> W2;
//				W2(0) = flow_data->W[(row-1)*3+0];
//				W2(1) = flow_data->W[(row-1)*3+1];
//				W2(2) = flow_data->W[(row-1)*3+2];
//				Matrix3<dreal> dWndW = inletBC_gradient(flo_opts, flow_data);
//				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, W2(0), W2(1), W2(2), W1(0), W1(1), W1(2));
//				d_lambdaw_dw *= flo_opts.scalar_d_eps;
//				Matrix3<dreal> dRdWn = -area_n * half * (analytical_flux_jacobian(gam, W2) - d_lambdaw_dw);
//				if(row==print_row) std::cout<<"dWndW "<<std::endl<<dWndW<<std::endl;
//				if(row==print_row) std::cout<<"dRdWn "<<std::endl<<dRdWn<<std::endl;
//				jacobian_diag = jacobian_diag + dRdWn*dWndW;
//			}
//			if (row == n_elem) { // Add boundary contribution
//				Vector3<dreal> W2;
//				W2(0) = flow_data->W[(row+1)*3+0];
//				W2(1) = flow_data->W[(row+1)*3+1];
//				W2(2) = flow_data->W[(row+1)*3+2];
//				Matrix3<dreal> dWpdW = outletBC_gradient(flo_opts, flow_data);
//				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, W1(0), W1(1), W1(2), W2(0), W2(1), W2(2));
//				d_lambdaw_dw *= flo_opts.scalar_d_eps;
//				Matrix3<dreal> dRdWp = area_p * half * (analytical_flux_jacobian(gam, W2) - d_lambdaw_dw);
//				//std::cout<<"Row "<<row<<std::endl<<jacobian_diag<<std::endl;
//				//std::cout<<"dWpdW "<<std::endl<<dWpdW<<std::endl;
//				//std::cout<<"dRdWp "<<std::endl<<dRdWp<<std::endl;
//				if(row==print_row) std::cout<<"dWpdW "<<std::endl<<dWpdW<<std::endl;
//				if(row==print_row) std::cout<<"dRdWp "<<std::endl<<dRdWp<<std::endl;
//				jacobian_diag = jacobian_diag + dRdWp*dWpdW;
//			}
//			if(row==print_row) std::cout<<"Row "<<row<<std::endl<<jacobian_diag<<std::endl;
//			if(print_jac) std::cout<<"Middle"<<std::endl;
//			if(print_jac) std::cout<<row<<std::endl<<jacobian_diag<<std::endl;
//			Matrix3<dreal> Middle;
//			for (int i = 0; i<3; i++) {
//				for (int j = 0; j<3; j++) {
//					Middle(i,j) = dRdW.coeff((row-1)*3+i,(row-1)*3+j);
//				}
//			}
//			if(print_jac) std::cout<<Middle<<std::endl<<std::endl;
//			if(print_jac) std::cout<<jacobian_diag-Middle<<std::endl<<std::endl;
//			if(print_jac) if((jacobian_diag-Middle).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;
//
//
//			// Add pseudo-timestep
//			dreal diag_identity = dx[row]/flow_data->dt[row];
//			jacobian_diag = jacobian_diag + diag_identity*identity;
//			if(row==print_row) std::cout<<"Row "<<row<<std::endl<<jacobian_diag<<std::endl;
//
//			if(row==print_row) {
//				std::vector<double> area_double(n_elem+1);
//				Flow_data<double> flow_data_double(n_elem);
//				for (int row = 0; row < n_elem+2; row++) {
//					for (int i_state = 0; i_state < 3; i_state++) {
//						flow_data_double.W[row*3+i_state] = valueof(flow_data->W[row*3+i_state]);
//					}
//				}
//				for (int row = 0; row < n_elem+1; row++) {
//					area_double[row] = valueof(area[row]);
//				}
//				Eigen::SparseMatrix<double> dRdW = eval_dRdW_dRdX_adolc(flo_opts, area_double, flow_data_double);
//				std::cout<<dRdW;
//				std::cout<<dRdW<<std::endl;
//				//for (int istate = 0; istate < 3; istate++) {
//					//for (int jstate = 0; istate < 3; istate++) {
//						//std::cout<<dRdW((row-1)*3+istate, (row-1)*3+jstate)<<" ";
//					//}
//					//std::cout<<std::endl;
//				//}
//			}
//
//			//jacobian_diag = 10*diag_identity*identity;
//
//
//			for (int i_state = 0; i_state < 3; i_state++) {
//				rhs(i_state) = -flow_data->residual[row*3+i_state];
//			}
//			if(row==print_row) std::cout<<"rhs "<<rhs<<std::endl<<std::endl;
//			for (int col = row-1; col < row; col++) {
//				Vector3<dreal> W2;
//				if (row!=n_elem) {
//					int i_p = row+1;
//					W2(0) = flow_data->W[i_p*3+0];
//					W2(1) = flow_data->W[i_p*3+1];
//					W2(2) = flow_data->W[i_p*3+2];
//					d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, W1(0), W1(1), W1(2), W2(0), W2(1), W2(2));
//					d_lambdaw_dw *= flo_opts.scalar_d_eps;
//					Matrix3<dreal> dRdWp = area_p * half * (analytical_flux_jacobian(gam, W2) - d_lambdaw_dw);
//
//					W2(0) = flow_data->W_stage1[i_p*3+0];
//					W2(1) = flow_data->W_stage1[i_p*3+1];
//					W2(2) = flow_data->W_stage1[i_p*3+2];
//
//					if(print_jac) std::cout<<"Upper"<<std::endl;
//					if(print_jac) std::cout<<dRdWp<<std::endl;
//					Matrix3<dreal> Upper;
//					for (int i = 0; i<3; i++) {
//						for (int j = 0; j<3; j++) {
//							Upper(i,j) = dRdW.coeff((row-1)*3+i,(i_p-1)*3+j);
//						}
//					}
//					if(print_jac) std::cout<<Upper<<std::endl<<std::endl;
//					if(print_jac) std::cout<<dRdWp-Upper<<std::endl<<std::endl;
//					if(print_jac) if((dRdWp-Upper).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;
//					if(print_jac) std::cout<<"W2 should be 0 on forward sweep"<<std::endl<<W2<<std::endl;
//
//
//					rhs -= dRdWp*W2;
//				}
//
//				if (col < 1) continue;
//
//				//const dreal normal = -1.0;
//				//W1(0) = flow_data->W[col*3+0] + flow_data->W_stage2[col*3+0];
//				//W1(1) = flow_data->W[col*3+1] + flow_data->W_stage2[col*3+1];
//				//W1(2) = flow_data->W[col*3+2] + flow_data->W_stage2[col*3+2];
//
//				//rhs = rhs - half*area[col]*normal*WtoF(gam,W1);
//
//				//W1(0) = flow_data->W[col*3+0];
//				//W1(1) = flow_data->W[col*3+1];
//				//W1(2) = flow_data->W[col*3+2];
//
//				//rhs = rhs + half*area[col]*normal*WtoF(gam,W1);
//
//				//const dreal u_j = W1(1)/W1(0);
//				//const dreal c_j = get_c(gam,W1(0),W1(1),W1(2));
//				//const dreal lambda = flo_opts.scalar_d_eps * (u_i+c_i + u_j+c_j)/2.0;
//
//				//for (int i_state = 0; i_state < 3; i_state++) {
//				//	rhs(i_state) = rhs(i_state) + half*lambda*flow_data->W_stage2[col*3+i_state]*area[col]*normal;
//				//}
//
//				W2(0) = flow_data->W[col*3+0];
//				W2(1) = flow_data->W[col*3+1];
//				W2(2) = flow_data->W[col*3+2];
//				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, W2(0), W2(1), W2(2), W1(0), W1(1), W1(2));
//				d_lambdaw_dw *= flo_opts.scalar_d_eps;
//				Matrix3<dreal> dRdWn = -area_n * half * (analytical_flux_jacobian(gam, W2) - d_lambdaw_dw);
//
//				W2(0) = flow_data->W_stage2[col*3+0];
//				W2(1) = flow_data->W_stage2[col*3+1];
//				W2(2) = flow_data->W_stage2[col*3+2];
//
//				rhs -= dRdWn*W2;
//
//				
//				if(print_jac) std::cout<<"Lower"<<std::endl;
//				if(print_jac) std::cout<<dRdWn<<std::endl;
//				Matrix3<dreal> Lower;
//				for (int i = 0; i<3; i++) {
//					for (int j = 0; j<3; j++) {
//						Lower(i,j) = dRdW.coeff((row-1)*3+i,(col-1)*3+j);
//					}
//				}
//				if(print_jac) std::cout<<Lower<<std::endl<<std::endl;
//				if(print_jac) std::cout<<dRdWn-Lower<<std::endl<<std::endl;
//				if(print_jac) if((dRdWn-Lower).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;
//				//abort();
//
//			}
//			Vector3<dreal> dW = jacobian_diag.fullPivLu().solve(rhs);
//			for (int i_state = 0; i_state < 3; i_state++) {
//				flow_data->W_stage2[row*3+i_state] = dW[i_state];
//			}
//			if(row==print_row) std::cout<<"dW* "<<dW<<std::endl<<std::endl;
//			if(row==print_row) std::cout<<"A "<<jacobian_diag<<std::endl<<std::endl;
//			if(row==print_row) std::cout<<"rhs "<<rhs<<std::endl<<std::endl;
//			//if(row==print_row) abort();
//		} // Row loop
//
//		// Backward sweep
//		dreal normdw = 0;
//		for (int row = n_elem; row > 0; row--) {
//
//			const dreal u_i = flow_data->W[row*3+1] / flow_data->W[row*3+0];
//			const dreal c_i = get_c(gam,flow_data->W[row*3+0],flow_data->W[row*3+1],flow_data->W[row*3+2]);
//
//			const int i_w_p = row+1;
//			const int i_w_n = row-1;
//
//			const dreal u_n = flow_data->W[i_w_n*3+1] / flow_data->W[i_w_n*3+0];
//			const dreal c_n = get_c(gam,flow_data->W[i_w_n*3+0],flow_data->W[i_w_n*3+1],flow_data->W[i_w_n*3+2]);
//			const dreal lambda_n = (u_n+c_n + u_i+c_i)/2.0;
//
//
//			const int i_face_p = row;
//			const int i_face_n = row-1;
//			const dreal area_p = area[i_face_p];
//			const dreal area_n = area[i_face_n];
//
//			//Vector3<dreal> W1;
//			//W1(0) = flow_data->W[row*3+0];
//			//W1(1) = flow_data->W[row*3+1];
//			//W1(2) = flow_data->W[row*3+2];
//			//Matrix3<dreal> jacobian_diag = analytical_flux_jacobian(gam, W1) * half;
//			//jacobian_diag(1,0) = jacobian_diag(1,0) + get_dpdw1(gam, W1(0), W1(1), W1(2));
//			//jacobian_diag(1,1) = jacobian_diag(1,1) + get_dpdw2(gam, W1(0), W1(1), W1(2));
//			//jacobian_diag(1,2) = jacobian_diag(1,2) + get_dpdw3(gam, W1(0), W1(1), W1(2));
//			//jacobian_diag = (area_p-area_n)*jacobian_diag;
//
//			//Matrix3<dreal> d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, flow_data->W[i_w_n*3+0], flow_data->W[i_w_n*3+1], flow_data->W[i_w_n*3+2], W1(0), W1(1), W1(2));
//			//jacobian_diag += half * flo_opts.scalar_d_eps * d_lambdaw_dw * area_n;
//			//d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, W1(0), W1(1), W1(2), flow_data->W[i_w_p*3+0], flow_data->W[i_w_p*3+1], flow_data->W[i_w_p*3+2]);
//			//jacobian_diag -= half * flo_opts.scalar_d_eps * d_lambdaw_dw * area_p;
//
//			Vector3<dreal> W1;
//			W1(0) = flow_data->W[row*3+0];
//			W1(1) = flow_data->W[row*3+1];
//			W1(2) = flow_data->W[row*3+2];
//			Matrix3<dreal> jacobian_diag = (analytical_flux_jacobian(gam, W1))*half;
//			jacobian_diag(1,0) = jacobian_diag(1,0) - get_dpdw1(gam, W1(0), W1(1), W1(2));
//			jacobian_diag(1,1) = jacobian_diag(1,1) - get_dpdw2(gam, W1(0), W1(1), W1(2));
//			jacobian_diag(1,2) = jacobian_diag(1,2) - get_dpdw3(gam, W1(0), W1(1), W1(2));
//			jacobian_diag = (area_p-area_n)*jacobian_diag;
//			
//			Matrix3<dreal> d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, flow_data->W[i_w_n*3+0], flow_data->W[i_w_n*3+1], flow_data->W[i_w_n*3+2], W1(0), W1(1), W1(2));
//			jacobian_diag += half * flo_opts.scalar_d_eps * d_lambdaw_dw * area_n;
//
//			d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, W1(0), W1(1), W1(2), flow_data->W[i_w_p*3+0], flow_data->W[i_w_p*3+1], flow_data->W[i_w_p*3+2]);
//			jacobian_diag -= half * flo_opts.scalar_d_eps * d_lambdaw_dw * area_p;
//
//
//			if (row == 1) { // Add boundary contribution
//				//Vector3<dreal> W2;
//				//W2(0) = flow_data->W[(row-1)*3+0];
//				//W2(1) = flow_data->W[(row-1)*3+1];
//				//W2(2) = flow_data->W[(row-1)*3+2];
//				//Matrix3<dreal> dWidWd = inletBC_gradient(flo_opts, flow_data);
//				//Matrix3<dreal> dRddWi = -half*(area_n) * (analytical_flux_jacobian(gam, W2) - lambda_n*identity);
//				//jacobian_diag = jacobian_diag + dRddWi*dWidWd;
//				Vector3<dreal> W2;
//				W2(0) = flow_data->W[(row-1)*3+0];
//				W2(1) = flow_data->W[(row-1)*3+1];
//				W2(2) = flow_data->W[(row-1)*3+2];
//				Matrix3<dreal> dWndW = inletBC_gradient(flo_opts, flow_data);
//				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, W2(0), W2(1), W2(2), W1(0), W1(1), W1(2));
//				d_lambdaw_dw *= flo_opts.scalar_d_eps;
//				Matrix3<dreal> dRdWn = -area_n * half * (analytical_flux_jacobian(gam, W2) - d_lambdaw_dw);
//				if(row==print_row) std::cout<<"dWndW "<<std::endl<<dWndW<<std::endl;
//				if(row==print_row) std::cout<<"dRdWn "<<std::endl<<dRdWn<<std::endl;
//				jacobian_diag = jacobian_diag + dRdWn*dWndW;
//			}
//			if (row == n_elem) { // Add boundary contribution
//				Vector3<dreal> W2;
//				W2(0) = flow_data->W[(row+1)*3+0];
//				W2(1) = flow_data->W[(row+1)*3+1];
//				W2(2) = flow_data->W[(row+1)*3+2];
//				Matrix3<dreal> dWodWd = outletBC_gradient(flo_opts, flow_data);
//				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, W1(0), W1(1), W1(2), W2(0), W2(1), W2(2));
//				d_lambdaw_dw *= flo_opts.scalar_d_eps;
//				Matrix3<dreal> dRdWp = area_p * half * (analytical_flux_jacobian(gam, W2) + d_lambdaw_dw);
//				jacobian_diag = jacobian_diag + dRdWp*dWodWd;
//			}
//
//			if(row==print_row) std::cout<<"Row "<<row<<std::endl<<jacobian_diag<<std::endl;
//			if(print_jac) std::cout<<"Middle"<<std::endl;
//			if(print_jac) std::cout<<row<<std::endl<<jacobian_diag<<std::endl;
//			Matrix3<dreal> Middle;
//			for (int i = 0; i<3; i++) {
//				for (int j = 0; j<3; j++) {
//					Middle(i,j) = dRdW.coeff((row-1)*3+i,(row-1)*3+j);
//				}
//			}
//			if(print_jac) std::cout<<Middle<<std::endl<<std::endl;
//			if(print_jac) std::cout<<jacobian_diag-Middle<<std::endl<<std::endl;
//			if(print_jac) if((jacobian_diag-Middle).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;
//
//			dreal diag_identity = dx[row]/flow_data->dt[row];
//			jacobian_diag = jacobian_diag + diag_identity*identity;
//
//			//jacobian_diag = 10*diag_identity*identity;
//
//			// Backward sweep
//			//rhs.setZero();
//			for (int i_state = 0; i_state < 3; i_state++) {
//				rhs(i_state) = -flow_data->residual[row*3+i_state];
//			}
//			for (int col = row+1; col < row+2; col++) {
//				Vector3<dreal> W2;
//				if (row!=1) {
//					int i_n = row-1;
//					W2(0) = flow_data->W[i_n*3+0];
//					W2(1) = flow_data->W[i_n*3+1];
//					W2(2) = flow_data->W[i_n*3+2];
//					d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 0, W2(0), W2(1), W2(2), W1(0), W1(1), W1(2));
//					d_lambdaw_dw *= flo_opts.scalar_d_eps;
//					Matrix3<dreal> dRdWn = -area_n * half * (analytical_flux_jacobian(gam, W2) - d_lambdaw_dw);
//
//					if(print_jac) std::cout<<"Lower"<<std::endl;
//					if(print_jac) std::cout<<dRdWn<<std::endl;
//					Matrix3<dreal> Lower;
//					for (int i = 0; i<3; i++) {
//						for (int j = 0; j<3; j++) {
//							Lower(i,j) = dRdW.coeff((row-1)*3+i,(i_n-1)*3+j);
//						}
//					}
//					if(print_jac) std::cout<<Lower<<std::endl<<std::endl;
//					if(print_jac) std::cout<<dRdWn-Lower<<std::endl<<std::endl;
//					if(print_jac) if((dRdWn-Lower).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;
//
//					W2(0) = flow_data->W_stage2[i_n*3+0];
//					W2(1) = flow_data->W_stage2[i_n*3+1];
//					W2(2) = flow_data->W_stage2[i_n*3+2];
//
//					rhs -= dRdWn*W2;
//				}
//				if (col > n_elem) continue;
//
//				//const dreal normal = 1.0;
//				//W1(0) = flow_data->W[col*3+0] + flow_data->W_stage1[col*3+0];
//				//W1(1) = flow_data->W[col*3+1] + flow_data->W_stage1[col*3+1];
//				//W1(2) = flow_data->W[col*3+2] + flow_data->W_stage1[col*3+2];
//
//				//rhs = rhs + half*area[col]*normal*WtoF(gam,W1);
//
//				//W1(0) = flow_data->W[col*3+0];
//				//W1(1) = flow_data->W[col*3+1];
//				//W1(2) = flow_data->W[col*3+2];
//
//				//rhs = rhs - half*area[col]*normal*WtoF(gam,W1);
//
//				//const dreal u_j = W1(1)/W1(0);
//				//const dreal c_j = get_c(gam,W1(0),W1(1),W1(2));
//				//const dreal lambda = flo_opts.scalar_d_eps * (u_i+c_i + u_j+c_j)/2.0;
//
//				//for (int i_state = 0; i_state < 3; i_state++) {
//				//	rhs(i_state) = rhs(i_state) + half*lambda*flow_data->W_stage1[col*3+i_state]*area[col]*normal;
//				//}
//				W2(0) = flow_data->W[col*3+0];
//				W2(1) = flow_data->W[col*3+1];
//				W2(2) = flow_data->W[col*3+2];
//				d_lambdaw_dw = get_d_lambdahalfw_dw(gam, 1, W1(0), W1(1), W1(2), W2(0), W2(1), W2(2));
//				d_lambdaw_dw *= flo_opts.scalar_d_eps;
//				Matrix3<dreal> dRdWp = area_p * half * (analytical_flux_jacobian(gam, W2) - d_lambdaw_dw);
//				if(print_jac) std::cout<<"Upper"<<std::endl;
//				if(print_jac) std::cout<<dRdWp<<std::endl;
//				Matrix3<dreal> Upper;
//				for (int i = 0; i<3; i++) {
//					for (int j = 0; j<3; j++) {
//						Upper(i,j) = dRdW.coeff((row-1)*3+i,(col-1)*3+j);
//					}
//				}
//				if(print_jac) std::cout<<Upper<<std::endl<<std::endl;
//				if(print_jac) std::cout<<dRdWp-Upper<<std::endl<<std::endl;
//				if(print_jac) if((dRdWp-Upper).norm() > 1e-13) std::cout<<"NOT OK****************"<<std::endl<<std::endl;
//
//				W2(0) = flow_data->W_stage1[col*3+0];
//				W2(1) = flow_data->W_stage1[col*3+1];
//				W2(2) = flow_data->W_stage1[col*3+2];
//				rhs -= dRdWp*W2;
//
//			}
//			Vector3<dreal> dW = jacobian_diag.fullPivLu().solve(rhs);
//			for (int i_state = 0; i_state < 3; i_state++) {
//				//flow_data->W_stage1[row*3+i_state] = flow_data->W_stage2[row*3+i_state] - dW[i_state];
//				flow_data->W_stage1[row*3+i_state] = dW[i_state];
//				normdw += pow(flow_data->W_stage1[row*3+i_state],2);
//			}
//		} // Row loop
//		std::cout<<"i_sweep"<<i_sweep<<" normdw"<<normdw<<std::endl;
//		if(print_jac) abort();
//	} // Sweep loop
//	for (int row = 1; row < n_elem+1; row++) {
//		for (int i_state = 0; i_state < 3; i_state++) {
//			flow_data->W[row*3+i_state] = flow_data->W[row*3+i_state] + flow_data->W_stage1[row*3+i_state];
//		}
//	}
//}
//template void lusgs( const struct Flow_options &flo_opts, const std::vector<double> &area, const std::vector<double> &dx, class Flow_data<double>* const flow_data);
////template void lusgs( const struct Flow_options &flo_opts, const std::vector<adouble> &area, const std::vector<double> &dx, class Flow_data<adouble>* const flow_data);
void lusgs(
	const struct Flow_options &flo_opts,
    const std::vector<adouble> &area,
    const std::vector<double> &dx,
    class Flow_data<adouble>* const flow_data)
{
	return;
}

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

