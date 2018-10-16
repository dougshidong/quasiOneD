// Time Stepping Schemes
#include "structures.h"
#include<vector>
#include<math.h>
#include<iostream>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include "flux.h"
#include "convert.h"
//#include "petscGMRES.h"
//#include "residuald1.h"
#include "boundary_conditions.h"

int ki, kip;

// Domain Residual R = FS_i+1/2 - FS_i-1/2 - Qi
void getDomainResi( 
	const struct Flow_options &flo_opts,
	const std::vector<double> &area,
	const std::vector<double> &W,
	std::vector<double>* const fluxes,
	std::vector<double>* const residual)
{
	getFlux(flo_opts, W, fluxes);
	const int n_elem = flo_opts.n_elem;
    for (int i_state = 0; i_state < 3; i_state++) {
        for (int i_res = 1; i_res < n_elem+1; i_res++) {
            const int kim = (i_res - 1) * 3 + i_state;
            const int ki  = (i_res + 0) * 3 + i_state;
            (*residual)[ki] = (*fluxes)[ki] * area[i_res] - (*fluxes)[kim] * area[i_res-1];

			if (i_state==1) { // Source Term
				double p = get_p(flo_opts.gam, W[i_res*3+0], W[i_res*3+1], W[i_res*3+2]);
				(*residual)[ki] -= p * (area[i_res] - area[i_res-1]);
			}
        }
    }
}

void EulerExplicitStep(
	const struct Flow_options &flo_opts,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    struct Flow_data* const flow_data);

void jamesonrk(
	const struct Flow_options &flo_opts,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    struct Flow_data* const flow_data);

//void eulerImplicit(
//    const std::vector<double> &area,
//    const std::vector<double> &dx,
//    const std::vector<double> &dt,
//    Flow_data &flow_data);

//void crankNicolson(
//    const std::vector<double> &area,
//    const std::vector<double> &dx,
//    const std::vector<double> &dt,
//    Flow_data &flow_data);

void stepInTime(
	const struct Flow_options &flo_opts,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    struct Flow_data* const flow_data)
{
	// Calculate Time Step
	int n_elem = flo_opts.n_elem;
	const int first_cell = 0;
	const int last_cell = n_elem+1;
	for (int i = 1; i < n_elem+1; i++) {
		const double u = flow_data->W[i*3+1] / flow_data->W[i*3+0];
		const double c = get_c(flo_opts.gam, flow_data->W[i*3+0], flow_data->W[i*3+1], flow_data->W[i*3+2]);
		flow_data->dt[i] = (flo_opts.CFL * dx[i]) / fabs(u + c);
	}
    flow_data->dt[first_cell] = flow_data->dt[first_cell+1];
    flow_data->dt[last_cell] = flow_data->dt[last_cell-1];

    if (flo_opts.time_scheme == 0) {
        EulerExplicitStep(flo_opts, area, dx, flow_data);
    } else if (flo_opts.time_scheme == 1) {
        abort();//rk4(area, dx, dt, flow_data);
    } else if (flo_opts.time_scheme == 2) {
        jamesonrk(flo_opts, area, dx, flow_data);
    } else if (flo_opts.time_scheme == 3) {
        abort();//eulerImplicit(area, dx, dt, flow_data);
    } else if (flo_opts.time_scheme == 4) {
        abort();//crankNicolson(area, dx, dt, flow_data);
    }
	inletBC(flo_opts, flow_data->dt[first_cell], dx[first_cell], flow_data);
	outletBC(flo_opts, flow_data->dt[last_cell], dx[last_cell], flow_data);
}

struct Flow_data stepInTime_noupdate(
	const struct Flow_options &flo_opts,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    const struct Flow_data &flow_data)
{
    struct Flow_data new_flow = flow_data;
	// Calculate Time Step
	int n_elem = flo_opts.n_elem;
	for (int i = 1; i < n_elem+1; i++) {
		double u = new_flow.W[i*3+1] / new_flow.W[i*3+0];
		double c = get_c(flo_opts.gam, new_flow.W[i*3+0], new_flow.W[i*3+1], new_flow.W[i*3+2]);
		new_flow.dt[i] = (flo_opts.CFL * dx[i]) / fabs(u + c);
	}
    if (flo_opts.time_scheme == 0) {
        EulerExplicitStep(flo_opts, area, dx, &new_flow);
    } else if (flo_opts.time_scheme == 1) {
        abort();//rk4(area, dx, dt, new_flow);
    } else if (flo_opts.time_scheme == 2) {
        jamesonrk(flo_opts, area, dx, &new_flow);
    } else if (flo_opts.time_scheme == 3) {
        abort();//eulerImplicit(area, dx, dt, new_flow);
    } else if (flo_opts.time_scheme == 4) {
        abort();//crankNicolson(area, dx, dt, new_flow);
    }

	const int first_cell = 0;
	const int last_cell = n_elem-1;
	inletBC(flo_opts, new_flow.dt[first_cell], dx[first_cell], &new_flow);
	outletBC(flo_opts, new_flow.dt[last_cell], dx[last_cell], &new_flow);

    return new_flow;
}


// Euler Explicit
void EulerExplicitStep(
	const struct Flow_options &flo_opts,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    struct Flow_data* const flow_data)
{
    getDomainResi(flo_opts, area, flow_data->W, &flow_data->fluxes, &flow_data->residual);

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
void jamesonrk(
	const struct Flow_options &flo_opts,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    struct Flow_data* const flow_data)
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
		getDomainResi(flo_opts, area, flow_data->W_stage, &flow_data->fluxes, &flow_data->residual);
        // Step in RK time
        for (int k = 0; k < 3; k++) {
            for (int i = 1; i < n_elem+1; i++) {
                ki = i * 3 + k;
                flow_data->W_stage[ki] = flow_data->W[ki] - (flow_data->dt[i-1] / (5 - r)) * flow_data->residual[ki] / dx[i-1];
            }
        }
    }

    for (int k = 0; k < 3; k++) {
        for (int i = 1; i < n_elem+1; i++) {
            ki = i * 3 + k;
            flow_data->residual[ki] = (flow_data->W_stage[ki] - flow_data->W[ki]) * dx[i] / flow_data->dt[i];
            flow_data->W[ki] = flow_data->W_stage[ki];
        }
    }
}

//void dFdW(
//    std::vector<double> &J,
//    double Sp, double Sm,
//    double rho, double u, double c)
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
//    double Sa = (Sp + Sm) / 2.0;
//    for (int i = 0; i < 9; i++)
//    {
//        J[i] *= Sa;
//    }
//    
//    double dpdw[3];
//    dpdw[0] = (gam - 1) / 2.0 * u * u;
//    dpdw[1] = - (gam - 1) * u;
//    dpdw[2] = (gam - 1);
//
//    J[3] -= dpdw[0] * (Sp - Sm);
//    J[4] -= dpdw[1] * (Sp - Sm);
//    J[5] -= dpdw[2] * (Sp - Sm);
//}

//void eulerImplicit(
//    const std::vector<double> &area,
//    const std::vector<double> &dx,
//    const std::vector<double> &dt,
//    Flow_data &flow_data)
//{
//    // Get flow_data.fluxes
//    getFlux(flow_data.fluxes, flow_data.W);
//    // Calculate Residuals
//    getDomainResi(flow_data, area);
//    Eigen::VectorXd RHS(3 * n_elem);
//    RHS.setZero();
//    Eigen::SparseMatrix<double> A(3 * n_elem, 3 * n_elem);
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
//    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
//    slusolver1.compute(A);
//    if (slusolver1.info() != 0)
//        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
//    Wt = slusolver1.solve(RHS);
//    double currentR = 0;
//    for (int i = 0; i < n_elem; i++)
//        currentR += flow_data.residual[i * 3 + 0] * flow_data.residual[i * 3 + 0];
//    currentR = sqrt(currentR);
//    double alpha = 1;
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
//    const std::vector<double> &area,
//    const std::vector<double> &dx,
//    const std::vector<double> &dt,
//    Flow_data &flow_data)
//{
//    // Get flow_data.fluxes
//    getFlux(flow_data.fluxes, flow_data.W);
//    // Calculate Residuals
//    getDomainResi(flow_data, area);
//    Eigen::VectorXd RHS(3 * n_elem);
//    RHS.setZero();
//    Eigen::SparseMatrix<double> A(3 * n_elem, 3 * n_elem);
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
//    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
//    slusolver1.compute(A);
//    if (slusolver1.info() != 0)
//        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
//    Wt = slusolver1.solve(RHS);
//    double currentR = 0;
//    for (int i = 0; i < n_elem; i++)
//        currentR += flow_data.residual[i * 3 + 0] * flow_data.residual[i * 3 + 0];
//    currentR = sqrt(currentR);
//    double alpha = 1;
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
//    const std::vector<double> &area,
//    const std::vector<double> &dx,
//    const std::vector<double> &dt,
//    Flow_data &flow_data);
//
// 4th order Runge - Kutta Stepping Scheme
// Need to figure out how to store more W
//void rk4(
//    const std::vector<double> &area,
//    const std::vector<double> &dx,
//    const std::vector<double> &dt,
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
//            Wtemp[ki] = ((double)1.0 / 6.0) * (flow_data.W[ki] + 2 * W1[ki] + 2 * W2[ki] + W3[ki]);
//            //flow_data.residual[ki] = (2 * Resi0[ki] + 2 * Resi1[ki] + Resi2[ki]) / 6.0;
//            flow_data.residual[ki] = (Wtemp[ki] - flow_data.W[ki]) * dx[i] / dt[i];
//            flow_data.W[ki] = Wtemp[ki];
//        }
//    }
//}

