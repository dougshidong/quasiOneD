#include<Eigen/Core>
#include<Eigen/Sparse>
#include<math.h>
#include<iostream>
#include"convert.h"
#include"flux.h"
#include"quasiOneD.h"
#include"flux.h"
#include"boundary_gradient.h"

using namespace Eigen;

void dFluxdW_scalard(
	const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector<double> &Ap_list,
    std::vector<double> &An_list);

std::vector<double> evaldQdW(
    std::vector<double> W,
    std::vector<double> area);

std::vector<double> evaldpdW(
	const double gamma,
    const std::vector<double> &W);

//************************************************************************
// Calculates Jacobian
SparseMatrix<double> evaldRdW(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data)
{
    // Do not use with implicit solvers
	if (flo_opts.time_scheme > 2) abort();
	// Only linearized scalar dissipation
    if (flo_opts.flux_scheme != 0) abort();

	int n_elem = flo_opts.n_elem;

    SparseMatrix<double> dRdW(3*n_elem, 3*n_elem);
    dRdW.reserve(27 * (n_elem - 2) + 27 * 4);

    // Get Jacobians and Fluxes
    std::vector<double> Ap(n_elem * 3 * 3, 0), An(n_elem * 3 * 3, 0);
	dFluxdW_scalard(flo_opts, flow_data.W, Ap, An);
	
    // Get Boundary Jacobians
    std::vector<double> dBidWi(9), dBidWd(9), dBodWd(9), dBodWo(9);
	dRdW_BC_inlet(flo_opts, flow_data.W, dBidWi, dBidWd);
	dRdW_BC_outlet(flo_opts, flow_data.W, dBodWd, dBodWo);

    // Evaluate dpdW
    std::vector<double> dpdW(3*n_elem, 0);
    dpdW = evaldpdW(flo_opts.gam, flow_data.W);

    // Input 4 lines where BC Jacobians occur
    // psi(1), psi(2), psi(n-1), psi(n)
    int Ri, Wi, k, rowi, coli;
    double val;
    for (int row = 0; row < 3; row++)
    {
        for (int col = 0; col < 3; col++)
        {
            k = row * 3 + col;
            // d(inlet)/d(inlet)
            // R0, W0
            Ri = 0; Wi = 0;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;

            val = - 1.0/n_elem/flow_data.dt[0] * dBidWi[k];
            dRdW.insert(rowi, coli) = val;

            // d(inlet)/d(domain)
            // R0, W1
            Ri = 0; Wi = 1;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;

            val = - 1.0/n_elem/flow_data.dt[0] *dBidWd[k];
            dRdW.insert(rowi, coli) = val;

            // d(outlet)/d(outlet)
            // R = n_elem - 1, W = n_elem - 1
            Ri = n_elem - 1; Wi = n_elem - 1;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;

            val = - 1e-0/n_elem/flow_data.dt[n_elem-1] *dBodWo[k];
            dRdW.insert(rowi, coli) = val;

            // d(outlet)/d(domain)
            // R = n_elem - 1, W = n_elem - 2
            Ri = n_elem - 1; Wi = n_elem - 2;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;

            val = - 1e-0/n_elem/flow_data.dt[n_elem-1] *dBodWd[k];
            dRdW.insert(rowi, coli) = val;
        }
    }
    for (int Ri = 1; Ri < n_elem - 1; Ri++)
    {
        Wi = Ri - 1;
        if (Wi >= 0)
        {
            for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++)
            {
                k = row * 3 + col;
                rowi = Ri * 3 + row;
                coli = Wi * 3 + col;

                val = - Ap[Wi * 9 + k] * area[Ri];
                dRdW.insert(rowi, coli) = val;
            }
        }

        Wi = Ri;
        if (Wi >= 0 && Wi <= n_elem - 1)
        {
            for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++)
            {
                k = row * 3 + col;
                rowi = Ri * 3 + row;
                coli = Wi * 3 + col;

                val = Ap[Wi * 9 + k] * area[Ri + 1];
                val -= An[Wi * 9 + k] * area[Ri];
                if (row == 1)
                {
                    val -= dpdW[Wi * 3 + col] * (area[Ri + 1] - area[Ri]);
                }

                dRdW.insert(rowi, coli) = val;
            }
        }

        Wi = Ri + 1;
        if (Wi <= n_elem - 1)
        {
            for (int row = 0; row < 3; row++)
            for (int col = 0; col < 3; col++)
            {
                k = row * 3 + col;
                rowi = Ri * 3 + row;
                coli = Wi * 3 + col;

                val = An[Wi * 9 + k] * area[Ri + 1];
                dRdW.insert(rowi, coli) = val;
            }
        }
    }
    return dRdW;
}

//SparseMatrix<double> evaldRdW_FD(
//    std::vector<double> W,
//    std::vector<double> area)
//{
//	int n_elem = W.size()/3;
//    SparseMatrix<double> dRdW(3 * n_elem, 3 * n_elem);
//    int Ri, Wi;
//    int rowi, coli;
//    std::vector<double> Wd(3 * n_elem, 0), Q(3 * n_elem, 0);
//    std::vector<double> Flux(3 * (n_elem + 1), 0);
//    std::vector<double> Resi1(3 * n_elem, 0), Resi2(3 * n_elem, 0);
//    std::vector<double> dRdW_block(9, 0), dRdWp(9, 0);
//    WtoQ(W, Q, area);
//    getFlux(Flux, W);
//    int ki, kip;
//    double pert;
//
//    // DR/DW
//    for (int Ri = 0; Ri < n_elem; Ri++) // LOOP OVER R
//    {
//        for (int Wi = 0; Wi < n_elem; Wi++) // LOOP OVER W
//        {
//            double h = 1e-5;
//            for (int statei = 0; statei < 3; statei++) // LOOP OVER STATEI
//            {
//                for (int i = 0; i < 3 * n_elem; i++)
//                    Wd[i] = W[i];
//
//                pert = W[Wi * 3 + statei] * h;
//                Wd[Wi * 3 + statei] = W[Wi * 3 + statei] + pert;
//
//                // RESI 1
//                // Inlet
//                if (Ri == 0) inletBC(Wd, Resi1, 1.0, 1.0);
//                // Outlet
//                else if (Ri == n_elem - 1) outletBC(Wd, Resi1, 1.0, 1.0);
//                // Domain
//                else
//                {
//                    WtoQ(Wd, Q, area);
//                    getFlux(Flux, Wd);
//
//                    for (int resii = 0; resii < 3; resii++)
//                    {
//                        ki = Ri * 3 + resii;
//                        kip = (Ri + 1) * 3 + resii;
//                        Resi1[ki] = Flux[kip] * area[Ri + 1] - Flux[ki] * area[Ri] - Q[ki];
//                    }
//                }
//
//                for (int i = 0; i < 3 * n_elem; i++)
//                    Wd[i] = W[i];
//
//                Wd[Wi * 3 + statei] = W[Wi * 3 + statei] - pert;
//                // RESI 2
//                // Inlet
//                if (Ri == 0) inletBC(Wd, Resi2, 1.0, 1.0);
//                // Outlet
//                else if (Ri == n_elem - 1) outletBC(Wd, Resi2, 1.0, 1.0);
//                // Domain
//                else
//                {
//                    WtoQ(Wd, Q, area);
//                    getFlux(Flux, Wd);
//
//                    for (int resii = 0; resii < 3; resii++)
//                    {
//                        ki = Ri * 3 + resii;
//                        kip = (Ri + 1) * 3 + resii;
//                        Resi2[ki] = Flux[kip] * area[Ri + 1] - Flux[ki] * area[Ri] - Q[ki];
//                    }
//                }
//
//                for (int resii = 0; resii < 3; resii++)
//                {
//                    ki = Ri * 3 + resii;
//                    dRdW_block[resii * 3 + statei] = (Resi1[ki] - Resi2[ki]) / (2 * pert);
//                }
//
//            } // END STATEI LOOP
//
//            for (int row = 0; row < 3; row++)
//            for (int col = 0; col < 3; col++)
//            {
//                rowi = Ri * 3 + row;
//                coli = Wi * 3 + col;
//                dRdW.insert(rowi, coli) = dRdW_block[row * 3 + col];
//            }
//        }  // END LOOP OVER W
//    } // END LOOP OVER R
//    return dRdW;
//}
//
//// Steger-Warming Flux Splitting
//void StegerJac(
//    std::vector<double> W,
//    std::vector<double> &Ap_list,
//    std::vector<double> &An_list,
//    std::vector<double> &Flux)
//{
//    double eps = 0.1;
//    double gam = 1.4;
//    double M[3][3] = {{0}},
//           Minv[3][3] = {{0}},
//           N[3][3] = {{0}},
//           Ninv[3][3] = {{0}},
//           lambdaP[3][3],
//           lambdaN[3][3];
//    double lambdaa[3];
//
//
//    double Ap[3][3], An[3][3], tempP[3][3], tempN[3][3], prefix[3][3], suffix[3][3];
//
//    std::vector<double> rho(n_elem), u(n_elem), p(n_elem), c(n_elem);
//    std::vector<double> Ap_list1(n_elem * 3 * 3, 0), An_list1(n_elem * 3 * 3, 0);
//
//    double beta = gam - 1;
//
//    for (int i = 0; i < n_elem; i++)
//    {
//        rho[i] = W[i * 3 + 0];
//        u[i] = W[i * 3 + 1] / rho[i];
//        p[i] = (gam-1) * (W[i * 3 + 2] - rho[i] * pow(u[i], 2) / 2);
//        c[i] = sqrt( gam * p[i] / rho[i] );
//    }
//
//
//    for (int i = 0; i < n_elem; i++)
//    {
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//        {
//            Ap[row][col] = 0;
//            An[row][col] = 0;
//            tempP[row][col] = 0;
//            tempN[row][col] = 0;
//            prefix[row][col] = 0;
//            suffix[row][col] = 0;
//            lambdaP[row][col] = 0;
//            lambdaN[row][col] = 0;
//        }
//
//        M[0][0] = 1.0;
//        M[1][0] = -u[i] / rho[i];
//        M[2][0] = 0.5 * u[i] * u[i] * beta;
//        M[1][1] = 1.0 / rho[i];
//        M[2][1] = -u[i] * beta;
//        M[2][2] = beta;
//        Minv[0][0] = 1.0;
//        Minv[1][0] = u[i];
//        Minv[2][0] = 0.5 * u[i] * u[i];
//        Minv[1][1] = rho[i];
//        Minv[2][1] = u[i] * rho[i];
//        Minv[2][2] = 1.0 / beta;
//        N[0][0] = 1.0;
//        N[1][1] = rho[i] * c[i];
//        N[2][1] = -rho[i] * c[i];
//        N[0][2] = -1.0 / (c[i] * c[i]);
//        N[1][2] = 1.0;
//        N[2][2] = 1.0;
//        Ninv[0][0] = 1.0;
//        Ninv[0][1] = 1.0 / (2.0 * c[i] * c[i]);
//        Ninv[0][2] = 1.0 / (2.0 * c[i] * c[i]);
//        Ninv[1][1] = 1.0 / (2.0 * rho[i] * c[i]);
//        Ninv[1][2] = -1.0 / (2.0 * rho[i] * c[i]);
//        Ninv[2][1] = 0.5;
//        Ninv[2][2] = 0.5;
//        lambdaa[0] = u[i];
//        lambdaa[1] = u[i] + c[i];
//        lambdaa[2] = u[i] - c[i];
//
//        for (int k = 0; k < 3; k++)
//            if (lambdaa[k] > 0)
//                lambdaP[k][k] = (lambdaa[k] + sqrt(pow(lambdaa[k], 2) + pow(eps, 2))) / 2.0;
//            else
//                lambdaN[k][k] = (lambdaa[k] - sqrt(pow(lambdaa[k], 2) + pow(eps, 2))) / 2.0;
//
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//            for (int k = 0; k < 3; k++)
//            {
//                prefix[row][col]+= Minv[row][k] * Ninv[k][col];
//                suffix[row][col]+= N[row][k] * M[k][col];
//            }
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//            for (int k = 0; k < 3; k++)
//            {
//                tempP[row][col] += prefix[row][k] * lambdaP[k][col];
//                tempN[row][col] += prefix[row][k] * lambdaN[k][col];
//            }
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//            for (int k = 0; k < 3; k++)
//            {
//                Ap[row][col]+= tempP[row][k] * suffix[k][col];
//                An[row][col]+= tempN[row][k] * suffix[k][col];
//            }
//        // could remove above loop and just use aplist and anlist
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//        {
//            int vec_pos = (i * 3 * 3) + (row * 3) + col;
//            Ap_list1[vec_pos] = Ap[row][col];
//            An_list1[vec_pos] = An[row][col];
//        }
//
//    }
//
//    for (int i = 1; i < n_elem; i++)
//    {
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//        {
//            int Ap_pos = ((i - 1) * 3 * 3) + (row * 3) + col;
//            int An_pos = (i * 3 * 3) + (row * 3) + col;
//            Flux[i * 3 + row] += Ap_list1[Ap_pos] * W[(i - 1) * 3 + col]
//                                 + An_list1[An_pos] * W[i * 3 + col];
//        }
//    }
//
//}

void JacobianCenter(
	const double gam,
    const double u,
	const double c,
    std::vector<double> &J)
{
    J[0] = 0.0;
    J[1] = 1.0;
    J[2] = 0.0;
    J[3] = u * u * (gam - 3.0) / 2.0;
    J[4] = u * (3.0 - gam);
    J[5] = gam - 1.0;
    J[6] = ( pow(u, 3) * (gam - 1.0) * (gam - 2.0) - 2.0 * u * c * c ) / (2.0 * (gam - 1.0));
    J[7] = ( 2.0 * c * c + u * u * ( -2.0 * gam * gam + 5.0 * gam - 3.0 ) )
           / (2.0 * (gam - 1.0));
    J[8] = u * gam;
}

std::vector<double> evaldlambdadW(
	const double gam,
    const double rho,
	const double rho_u,
	const double e)
{
    std::vector<double> dlambdadWp(3), dlambdadW(3);
    const double p = get_p(gam, rho, rho_u, e);//(gam - 1) * ( e - rho_u*rho_u/rho / 2 );
    const double c = get_c(gam, rho, rho_u, e);//sqrt( gam * p / rho );

    // dlambda/dWp
    const double dlambdadr = -c / (4.0 * rho);
    const double dlambdadu = 0.5;
    const double dlambdadp = c / (4.0 * p);

    dlambdadWp[0] = dlambdadr;
    dlambdadWp[1] = dlambdadu;
    dlambdadWp[2] = dlambdadp;

    std::vector<double> dWpdW(9, 0);
    eval_dWpdW(gam, rho, rho_u, &dWpdW);
    dlambdadW[0] = 0;
    dlambdadW[1] = 0;
    dlambdadW[2] = 0;
    for (int row = 0; row < 1; row++)
    for (int col = 0; col < 3; col++)
    for (int k = 0; k < 3; k++)
        dlambdadW[row * 3 + col] += dlambdadWp[row * 3 + k] * dWpdW[k * 3 + col];

    return dlambdadW;
}

void dFluxdW_scalard(
	const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector<double> &Ap_list,
    std::vector<double> &An_list)
{
	const int n_elem = flo_opts.n_elem;
	const double gam = flo_opts.gam;
	const double scalar_d_eps = flo_opts.scalar_d_eps;

    std::vector<double> u(n_elem), c(n_elem);
	for (int i=0; i<n_elem; i++) {
		u[i] = W[i*3+1] / W[i*3+0];
		c[i] = get_c(gam, W[i*3+0], W[i*3+1], W[i*3+2]);
	}

    int vec_pos, k;
    double lamb;

    std::vector<double> J(9, 0);
    std::vector<double> dlambdaPdW(3, 0);
    std::vector<double> dlambdaNdW(3, 0);
    // A+
    for (int i = 0; i < n_elem - 1; i++)
    {
        // dF/dW
        JacobianCenter(gam, u[i], c[i], J);

        // lambda
        lamb = (u[i] + u[i + 1] + c[i] + c[i + 1]) / 2.0;
//      // dlambdaP/dW
        dlambdaPdW = evaldlambdadW(gam, W[i*3+0], W[i*3+1], W[i*3+2]);

        for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++)
        {
            vec_pos = (i * 9) + (row * 3) + col; // NOT Transposed
            k = row * 3 + col;
            Ap_list[vec_pos] = J[k] / 2.0 - dlambdaPdW[col] * scalar_d_eps
                               * (W[(i + 1) * 3 + row] - W[i * 3 + row]) / 2.0;
            if (row == col)
            {
                Ap_list[vec_pos] += scalar_d_eps * lamb / 2.0;
            }
        }
    }

    // A-
    for (int i = 1; i < n_elem; i++)
    {
        // dF/dW
        JacobianCenter(gam, u[i], c[i], J);

        // lambda
        lamb = (u[i] + u[i - 1] + c[i] + c[i - 1]) / 2.0;
//      // dlambdaN/dW
        dlambdaNdW = evaldlambdadW(gam, W[i*3+0], W[i*3+1], W[i*3+2]);

        for (int row = 0; row < 3; row++)
        for (int col = 0; col < 3; col++)
        {
            vec_pos = (i * 9) + (row * 3) + col; // NOT Transposed
            k = row * 3 + col;
            An_list[vec_pos] = J[k] / 2.0 - dlambdaNdW[col] * scalar_d_eps
                               * (W[i * 3 + row] - W[(i - 1) * 3 + row]) / 2.0;
            if (row == col)
            {
                An_list[vec_pos] -= scalar_d_eps * lamb / 2.0;
            }
        }
    }
}

//std::vector<double> evaldQdW(
//    std::vector<double> W,
//    std::vector<double> area)
//{
//    double dpdw[3], rho, u, dArea;
//	int n_elem = W.size()/3;
//    std::vector<double> dQdW(3 * n_elem);
//    for (int i = 0; i < n_elem; i++)
//    {
//        rho = W[i * 3 + 0];
//        u = W[i * 3 + 1] / rho;
//
//        dpdw[0] = (gam - 1) / 2.0 * u * u;
//        dpdw[1] = - (gam - 1) * u;
//        dpdw[2] = (gam - 1);
//
//        dArea = area[i + 1] - area[i];
//
//        dQdW[i * 3 + 0] = dpdw[0] * dArea;
//        dQdW[i * 3 + 1] = dpdw[1] * dArea;
//        dQdW[i * 3 + 2] = dpdw[2] * dArea;
//    }
//    return dQdW;
//}

std::vector<double> evaldpdW(
	const double gam,
    const std::vector<double> &W)
{
	const int n_elem = W.size()/3;
    std::vector<double> dpdW(3 * n_elem);
    for (int i = 0; i < n_elem; i++)
    {
        double rho = W[i * 3 + 0];
        double u = W[i * 3 + 1] / rho;

        dpdW[i * 3 + 0] = (gam - 1) / 2.0 * u * u;
        dpdW[i * 3 + 1] = - (gam - 1) * u;
        dpdW[i * 3 + 2] = (gam - 1);
    }
    return dpdW;
}
// DERIVATIVES WRT GEOMETRY

MatrixXd evaldRdArea(
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data)
{
	const int n_elem = flo_opts.n_elem;
	std::vector<double> fluxes(3*(n_elem+1));
	getFlux(flo_opts, flow_data.W, fluxes);

    MatrixXd dRdArea(3 * n_elem, n_elem + 1);

    std::vector<double> Q(3 * n_elem, 0), p(n_elem);
    int Si, kR, kS;
    dRdArea.setZero();
    for (int Ri = 1; Ri < n_elem - 1; Ri++)
    {
		const double p = get_p(flo_opts.gam, flow_data.W[Ri*3+0], flow_data.W[Ri*3+1], flow_data.W[Ri*3+2]);
        for (int k = 0; k < 3; k++)
        {
            kR = Ri * 3 + k;

            Si = Ri;
            kS = Si * 3 + k;
            dRdArea(kR, Si) = -fluxes[kS];
            if (k == 1) dRdArea(kR, Si) += p;

            Si = Ri + 1;
            kS = Si * 3 + k;
            dRdArea(kR, Si) = fluxes[kS];
            if (k == 1) dRdArea(kR, Si) += -p;
        }
    }
    return dRdArea;
}

//MatrixXd evaldRdArea_FD(
//    std::vector<double> Flux,
//    std::vector<double> area,
//    std::vector<double> W)
//{
//    MatrixXd dRdArea(3 * n_elem, n_elem + 1);
//    dRdArea.setZero();
//    std::vector<double> Resi0(3 * n_elem, 0), Resi1(3 * n_elem, 0), Resi2(3 * n_elem, 0);
//    std::vector<double> Sd(n_elem + 1, 0);
//    std::vector<double> Q(3 * n_elem, 0);
//    double h = 0.000000001;
//    double pert;
//    int ki, kip;
//    for (int Ri = 1; Ri < n_elem - 1; Ri++)
//    {
//        for (int Si = 0; Si < n_elem + 1; Si++)
//        {
//            for (int m = 0; m < n_elem + 1; m++)
//                Sd[m] = area[m];
//
//            pert = area[Si] * h;
//            Sd[Si] = area[Si] + pert;
//
//            WtoQ(W, Q, Sd);
//
//            for (int k = 0; k < 3; k++)
//            {
//                ki = Ri * 3 + k;
//                kip = (Ri + 1) * 3 + k;
//                Resi1[ki] = Flux[kip] * Sd[Ri + 1] - Flux[ki] * Sd[Ri] - Q[ki];
//            }
//
//            for (int m = 0; m < n_elem + 1; m++)
//                Sd[m] = area[m];
//
//            Sd[Si] = area[Si] - pert;
//
//            WtoQ(W, Q, Sd);
//
//            for (int k = 0; k < 3; k++)
//            {
//                ki = Ri * 3 + k;
//                kip = (Ri + 1) * 3 + k;
//                Resi2[ki] = Flux[kip] * Sd[Ri + 1] - Flux[ki] * Sd[Ri] - Q[ki];
//                dRdArea(ki, Si) = (Resi1[ki] - Resi2[ki]) / (2 * pert);
//            }
//        }
//    }
//
//    return dRdArea;
//}
