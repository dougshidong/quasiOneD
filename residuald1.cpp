#include<Eigen/Core>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<math.h>
#include<iostream>
#include"convert.h"
#include"flux.h"
#include"quasiOneD.h"
#include"flux.h"
#include"boundary_gradient.h"
#include"boundary_conditions.h"
#include"timestep.h"


using namespace Eigen;

void dFluxdW_scalard(
	const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector<double> &dFluxdWL,
    std::vector<double> &dFluxdWR);

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

	const int n_elem = flo_opts.n_elem;
    const int n_resi = n_elem*3;
	const int n_face = n_elem+1;
    const int n_flux = n_face*3;

    SparseMatrix<double> dRdW(n_resi, n_resi);
    dRdW.reserve(9*n_resi);

    //const int n_rows = 3*(n_elem-2);
    //const int n_stencil = 3;
    //const int n_state = 3;
    //SparseMatrix<double> dRdW(n_rows, n_rows);
    //dRdW.reserve(n_stencil*n_state*n_rows);

    // Get Jacobians and Fluxes
    std::vector<double> dFluxdWLeft(n_face * 3 * 3, 0), dFluxdWRight(n_face * 3 * 3, 0);
	dFluxdW_scalard(flo_opts, flow_data.W, dFluxdWLeft, dFluxdWRight);
	
    // Evaluate dpdW
    std::vector<double> dpdW = evaldpdW(flo_opts.gam, flow_data.W);

    // Input 4 lines where BC Jacobians occur
    // psi(1), psi(2), psi(n-1), psi(n)
    int Ri, Wi;
    for (int Ri = 1; Ri < n_elem + 1; Ri++) {
        // dR/dWLeft
        const int iFace_left = Ri-1;
        const int iFace_right = Ri;

        Wi = Ri - 1;
        if (Wi > 0) {
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    const int k = row * 3 + col;

                    const double val = - dFluxdWLeft.at(iFace_left * 9 + k) * area.at(iFace_left);

                    const int irow_glob = (Ri-1) * 3 + row;
                    const int icol_glob = (Wi-1) * 3 + col;
                    dRdW.insert(irow_glob, icol_glob) = val;
                }
            }
        }

        // dR/dWCenter
        Wi = Ri;
        if (Wi > 0 && Wi < n_elem + 1) {
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    const int k = row * 3 + col;

                    double val = dFluxdWLeft.at(iFace_right * 9 + k) * area.at(iFace_right);
                    val -= dFluxdWRight.at(iFace_left * 9 + k) * area.at(iFace_left);
                    if (row == 1) val -= dpdW.at(Wi * 3 + col) * (area.at(iFace_right) - area.at(iFace_left));

                    const int irow_glob = (Ri-1) * 3 + row;
                    const int icol_glob = (Wi-1) * 3 + col;
                    dRdW.insert(irow_glob, icol_glob) = val;
                }
            }
        }

        // dR/dWRight
        Wi = Ri + 1;
        if (Wi < n_elem + 1) {
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    const int k = row * 3 + col;
                    const double val = dFluxdWRight.at(iFace_right * 9 + k) * area.at(iFace_right);

                    const int irow_glob = (Ri-1) * 3 + row;
                    const int icol_glob = (Wi-1) * 3 + col;
                    dRdW.insert(irow_glob, icol_glob) = val;
                }
            }
        }
    }

    // Get Boundary Jacobians
    std::vector<double> dBidWi(9), dBidWd(9), dBodWd(9), dBodWo(9);
    std::vector<double> dFluxdW(9), dFluxdW_dWdW(9);
	dRdW_BC_inlet(flo_opts, flow_data.W, dBidWi, dBidWd);
	dRdW_BC_outlet(flo_opts, flow_data.W, dBodWd, dBodWo);

    MatrixXd dFluxdW_e(3,3), dBidWi_e(3,3), dBidWd_e(3,3), dBodWd_e(3,3), dBodWo_e(3,3), dFluxdW_dWdW_e(3,3);
    const int Ri_first = 1;
    const int Wi_first = 1;
    const int iFace_left = Ri_first-1;
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            const int k = row * 3 + col;
            dFluxdW.at(k) = -dFluxdWLeft.at(iFace_left * 9 + k) * area.at(iFace_left);
            dFluxdW_dWdW.at(k) = 0;
            dBidWi_e(row,col) = dBidWi.at(k);
            dBidWd_e(row,col) = dBidWd.at(k);
            dBodWo_e(row,col) = dBodWo.at(k);
            dBodWd_e(row,col) = dBodWd.at(k);
            dFluxdW_e(row, col) = dFluxdW.at(k);
        }
    }

    dFluxdW_dWdW_e.setZero();

    MatrixXd identity = MatrixXd::Identity(3,3);
    MatrixXd actual_dWinlet_dWdomain(3,3);
    if ((identity-dBidWi_e).norm() == 0) {
        dFluxdW_dWdW_e.setZero();
    } else {
        MatrixXd actual_dWinlet_dWdomain = ((identity - dBidWi_e).inverse())*dBidWd_e;
    }
    dFluxdW_dWdW_e = dFluxdW_e * actual_dWinlet_dWdomain;

    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            const int irow_glob = (Ri_first-1) * 3 + row;
            const int icol_glob = (Wi_first-1) * 3 + col;
            const int k = row * 3 + col;
            dRdW.coeffRef(irow_glob, icol_glob) += dFluxdW_dWdW_e(row, col);
        }
    }

    const int Ri_last = n_elem;
    const int Wi_last = Ri_last;
    const int iFace_right = Ri_last-1;
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            const int k = row * 3 + col;
            dFluxdW.at(k) = dFluxdWRight.at(iFace_right * 9 + k) * area.at(iFace_right);
            dFluxdW_dWdW.at(k) = 0;
            dFluxdW_e(row, col) = dFluxdW.at(k);
        }
    }
    for (int row = 0; row < 3; row++)
    for (int col = 0; col < 3; col++)
    for (int k = 0; k < 3; k++)
        dFluxdW_dWdW.at(row * 3 + col) += dFluxdW.at(row * 3 + k) * dBodWd.at(k * 3 + col);

    dFluxdW_dWdW_e.setZero();
    dFluxdW_dWdW_e = dFluxdW_e * (dBodWd_e + dBodWo_e*dBodWd_e);
    MatrixXd actual_dWoutlet_dWdomain = ((identity - dBodWo_e).inverse())*dBodWd_e;
    dFluxdW_dWdW_e = dFluxdW_e * actual_dWoutlet_dWdomain;
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            const int irow_glob = (Ri_last-1) * 3 + row;
            const int icol_glob = (Wi_last-1) * 3 + col;
            const int k = row * 3 + col;
            dRdW.coeffRef(irow_glob, icol_glob) += dFluxdW_dWdW_e(row, col);
        }
    }

    return dRdW;
}

//
//// Steger-Warming Flux Splitting
//void StegerJac(
//    std::vector<double> W,
//    std::vector<double> &dFluxdWL,
//    std::vector<double> &dFluxdWR,
//    std::vector<double> &Flux)
//{
//    double eps = 0.1;
//    double gam = 1.4;
//    double M.at(3).at(3) = {{0}},
//           Minv.at(3).at(3) = {{0}},
//           N.at(3).at(3) = {{0}},
//           Ninv.at(3).at(3) = {{0}},
//           lambdaP.at(3).at(3),
//           lambdaN.at(3).at(3);
//    double lambdaa.at(3);
//
//
//    double Ap.at(3).at(3), An.at(3).at(3), tempP.at(3).at(3), tempN.at(3).at(3), prefix.at(3).at(3), suffix.at(3).at(3);
//
//    std::vector<double> rho(n_elem), u(n_elem), p(n_elem), c(n_elem);
//    std::vector<double> Ap_list1(n_elem * 3 * 3, 0), An_list1(n_elem * 3 * 3, 0);
//
//    double beta = gam - 1;
//
//    for (int i = 0; i < n_elem; i++)
//    {
//        rho.at(i) = W.at(i * 3 + 0);
//        u.at(i) = W.at(i * 3 + 1) / rho.at(i);
//        p.at(i) = (gam-1) * (W.at(i * 3 + 2) - rho.at(i) * pow(u.at(i), 2) / 2);
//        c.at(i) = sqrt( gam * p.at(i) / rho.at(i) );
//    }
//
//
//    for (int i = 0; i < n_elem; i++)
//    {
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//        {
//            Ap.at(row).at(col) = 0;
//            An.at(row).at(col) = 0;
//            tempP.at(row).at(col) = 0;
//            tempN.at(row).at(col) = 0;
//            prefix.at(row).at(col) = 0;
//            suffix.at(row).at(col) = 0;
//            lambdaP.at(row).at(col) = 0;
//            lambdaN.at(row).at(col) = 0;
//        }
//
//        M.at(0).at(0) = 1.0;
//        M.at(1).at(0) = -u.at(i) / rho.at(i);
//        M.at(2).at(0) = 0.5 * u.at(i) * u.at(i) * beta;
//        M.at(1).at(1) = 1.0 / rho.at(i);
//        M.at(2).at(1) = -u.at(i) * beta;
//        M.at(2).at(2) = beta;
//        Minv.at(0).at(0) = 1.0;
//        Minv.at(1).at(0) = u.at(i);
//        Minv.at(2).at(0) = 0.5 * u.at(i) * u.at(i);
//        Minv.at(1).at(1) = rho.at(i);
//        Minv.at(2).at(1) = u.at(i) * rho.at(i);
//        Minv.at(2).at(2) = 1.0 / beta;
//        N.at(0).at(0) = 1.0;
//        N.at(1).at(1) = rho.at(i) * c.at(i);
//        N.at(2).at(1) = -rho.at(i) * c.at(i);
//        N.at(0).at(2) = -1.0 / (c.at(i) * c.at(i));
//        N.at(1).at(2) = 1.0;
//        N.at(2).at(2) = 1.0;
//        Ninv.at(0).at(0) = 1.0;
//        Ninv.at(0).at(1) = 1.0 / (2.0 * c.at(i) * c.at(i));
//        Ninv.at(0).at(2) = 1.0 / (2.0 * c.at(i) * c.at(i));
//        Ninv.at(1).at(1) = 1.0 / (2.0 * rho.at(i) * c.at(i));
//        Ninv.at(1).at(2) = -1.0 / (2.0 * rho.at(i) * c.at(i));
//        Ninv.at(2).at(1) = 0.5;
//        Ninv.at(2).at(2) = 0.5;
//        lambdaa.at(0) = u.at(i);
//        lambdaa.at(1) = u.at(i) + c.at(i);
//        lambdaa.at(2) = u.at(i) - c.at(i);
//
//        for (int k = 0; k < 3; k++)
//            if (lambdaa.at(k) > 0)
//                lambdaP.at(k).at(k) = (lambdaa.at(k) + sqrt(pow(lambdaa.at(k), 2) + pow(eps, 2))) / 2.0;
//            else
//                lambdaN.at(k).at(k) = (lambdaa.at(k) - sqrt(pow(lambdaa.at(k), 2) + pow(eps, 2))) / 2.0;
//
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//            for (int k = 0; k < 3; k++)
//            {
//                prefix.at(row).at(col)+= Minv.at(row).at(k) * Ninv.at(k).at(col);
//                suffix.at(row).at(col)+= N.at(row).at(k) * M.at(k).at(col);
//            }
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//            for (int k = 0; k < 3; k++)
//            {
//                tempP.at(row).at(col) += prefix.at(row).at(k) * lambdaP.at(k).at(col);
//                tempN.at(row).at(col) += prefix.at(row).at(k) * lambdaN.at(k).at(col);
//            }
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//            for (int k = 0; k < 3; k++)
//            {
//                Ap.at(row).at(col)+= tempP.at(row).at(k) * suffix.at(k).at(col);
//                An.at(row).at(col)+= tempN.at(row).at(k) * suffix.at(k).at(col);
//            }
//        // could remove above loop and just use aplist and anlist
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//        {
//            int vec_pos = (i * 3 * 3) + (row * 3) + col;
//            Ap_list1.at(vec_pos) = Ap.at(row).at(col);
//            An_list1.at(vec_pos) = An.at(row).at(col);
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
//            Flux.at(i * 3 + row) += Ap_list1.at(Ap_pos) * W.at((i - 1) * 3 + col)
//                                 + An_list1.at(An_pos) * W.at(i * 3 + col);
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
    J.at(0) = 0.0;
    J.at(1) = 1.0;
    J.at(2) = 0.0;
    J.at(3) = u * u * (gam - 3.0) / 2.0;
    J.at(4) = u * (3.0 - gam);
    J.at(5) = gam - 1.0;
    J.at(6) = ( pow(u, 3) * (gam - 1.0) * (gam - 2.0) - 2.0 * u * c * c ) / (2.0 * (gam - 1.0));
    J.at(7) = ( 2.0 * c * c + u * u * ( -2.0 * gam * gam + 5.0 * gam - 3.0 ) )
           / (2.0 * (gam - 1.0));
    J.at(8) = u * gam;
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

    dlambdadWp.at(0) = dlambdadr;
    dlambdadWp.at(1) = dlambdadu;
    dlambdadWp.at(2) = dlambdadp;

    std::vector<double> dWpdW(9, 0);
    eval_dWpdW(gam, rho, rho_u, &dWpdW);
    dlambdadW.at(0) = 0;
    dlambdadW.at(1) = 0;
    dlambdadW.at(2) = 0;
    for (int row = 0; row < 1; row++)
    for (int col = 0; col < 3; col++)
    for (int k = 0; k < 3; k++)
        dlambdadW.at(row * 3 + col) += dlambdadWp.at(row * 3 + k) * dWpdW.at(k * 3 + col);

    return dlambdadW;
}

void dFluxdW_scalard(
	const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector<double> &dFluxdWL,
    std::vector<double> &dFluxdWR)
{
	const int n_elem = flo_opts.n_elem;
    const int n_resi = n_elem*3;
	const int n_face = n_elem+1;
    const int n_flux = n_face*3;

	const double gam = flo_opts.gam;
	const double scalar_d_eps = flo_opts.scalar_d_eps;

    std::vector<double> u(n_elem+2), c(n_elem+2);
	for (int i=0; i<n_elem+2; i++) {
		u.at(i) = W.at(i*3+1) / W.at(i*3+0);
		c.at(i) = get_c(gam, W.at(i*3+0), W.at(i*3+1), W.at(i*3+2));
	}

    int vec_pos, k;
    double lamb;

    std::vector<double> J(9, 0);
    std::vector<double> dLambdadW(3, 0);

    // dFlux/dWLeft
    // w0    |      w1   |   ...   |     WL    |    WR     |
    //    i_face=0              i_face-1    i_face      i_face+1
    for (int i_face = 0; i_face < n_face; i_face++) {
        const int iw_left = i_face;
        const int iw_right = i_face+1;
        // lambda
        lamb = (u.at(iw_left) + u.at(iw_right) + c.at(iw_left) + c.at(iw_right)) / 2.0;

        // dFLeft/dWLeft
        JacobianCenter(gam, u.at(iw_left), c.at(iw_left), J);
        // dlambda/dWLeft
        dLambdadW = evaldlambdadW(gam, W.at(iw_left*3+0), W.at(iw_left*3+1), W.at(iw_left*3+2));

        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                vec_pos = (i_face * 9) + (row * 3) + col; // NOT Transposed
                k = row * 3 + col;
                dFluxdWL.at(vec_pos) = J.at(k) / 2.0 - dLambdadW.at(col) * scalar_d_eps * (W.at(iw_right * 3 + row) - W.at(iw_left * 3 + row)) / 2.0;
                if (row == col) dFluxdWL.at(vec_pos) += scalar_d_eps * lamb / 2.0;
            }
        }

        // dFLeft/dWRight
        JacobianCenter(gam, u.at(iw_right), c.at(iw_right), J);
        // dlambda/dWRight
        dLambdadW = evaldlambdadW(gam, W.at(iw_right*3+0), W.at(iw_right*3+1), W.at(iw_right*3+2));
        for (int row = 0; row < 3; row++) { 
            for (int col = 0; col < 3; col++) {
                vec_pos = (i_face * 9) + (row * 3) + col; // NOT Transposed
                k = row * 3 + col;
                dFluxdWR.at(vec_pos) = J.at(k) / 2.0 - dLambdadW.at(col) * scalar_d_eps * (W.at(iw_right * 3 + row) - W.at(iw_left * 3 + row)) / 2.0;
                if (row == col) dFluxdWR.at(vec_pos) -= scalar_d_eps * lamb / 2.0;
            }
        }
    }

}

//std::vector<double> evaldQdW(
//    std::vector<double> W,
//    std::vector<double> area)
//{
//    double dpdw.at(3), rho, u, dArea;
//	int n_elem = W.size()/3;
//    std::vector<double> dQdW(n_resi);
//    for (int i = 0; i < n_elem; i++)
//    {
//        rho = W.at(i * 3 + 0);
//        u = W.at(i * 3 + 1) / rho;
//
//        dpdw.at(0) = (gam - 1) / 2.0 * u * u;
//        dpdw.at(1) = - (gam - 1) * u;
//        dpdw.at(2) = (gam - 1);
//
//        dArea = area.at(i + 1) - area.at(i);
//
//        dQdW.at(i * 3 + 0) = dpdw.at(0) * dArea;
//        dQdW.at(i * 3 + 1) = dpdw.at(1) * dArea;
//        dQdW.at(i * 3 + 2) = dpdw.at(2) * dArea;
//    }
//    return dQdW;
//}

std::vector<double> evaldpdW(
	const double gam,
    const std::vector<double> &W)
{
	const int n_w = W.size()/3;
    std::vector<double> dpdW(W.size());
    for (int i = 0; i < n_w; i++)
    {
        double rho = W.at(i * 3 + 0);
        double u = W.at(i * 3 + 1) / rho;

        dpdW.at(i * 3 + 0) = (gam - 1) / 2.0 * u * u;
        dpdW.at(i * 3 + 1) = - (gam - 1) * u;
        dpdW.at(i * 3 + 2) = (gam - 1);
    }
    return dpdW;
}
// DERIVATIVES WRT GEOMETRY

MatrixXd evaldRdArea(
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data)
{
	const int n_elem = flo_opts.n_elem;
    const int n_resi = n_elem*3;
	const int n_face = n_elem+1;
	const int n_flux = n_face*3;
	std::vector<double> fluxes(n_flux, 0);
	getFlux(flo_opts, flow_data.W, &fluxes);

    MatrixXd dRdArea(n_resi, n_face);

    std::vector<double> Q(n_resi, 0);
    int Si, kS;
    dRdArea.setZero();
    for (int Ri = 1; Ri < n_elem+1; Ri++)
    {
		const double p = get_p(flo_opts.gam, flow_data.W.at(Ri*3+0), flow_data.W.at(Ri*3+1), flow_data.W.at(Ri*3+2));
        for (int k = 0; k < 3; k++) {
            const int irow_glob = (Ri-1) * 3 + k;

            int i_area = Ri-1;
            int i_flux = i_area * 3 + k;
            dRdArea(irow_glob, i_area) = -fluxes.at(i_flux);
            if (k == 1) dRdArea(irow_glob, i_area) += p;

            i_area = Ri;
            i_flux = i_area * 3 + k;
            dRdArea(irow_glob, i_area) = fluxes.at(i_flux);
            if (k == 1) dRdArea(irow_glob, i_area) += -p;
        }
    }
    return dRdArea;
}

MatrixXd evaldRdArea_FD(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data)
{
	const int n_elem = flo_opts.n_elem;
    const int n_resi = n_elem*3;
	const int n_face = n_elem+1;
    MatrixXd dRdArea(n_resi, n_face);
    dRdArea.setZero();

    std::vector<double> pert_area1 = area;
    std::vector<double> pert_area2 = area;
    struct Flow_data pert_flow1 = flow_data;
    struct Flow_data pert_flow2 = flow_data;

    double dh = 1e-03;
    for (int Ri = 1; Ri < n_elem - 1; Ri++) {
        for (int Si = 0; Si < n_elem + 1; Si++) {
            const double pertA = area.at(Si) * dh;

            pert_area1 = area;
            pert_flow1 = flow_data;
            pert_area1.at(Si) = area.at(Si) + pertA;

            getDomainResi(flo_opts, pert_area1, pert_flow1.W, &(pert_flow1.fluxes), &(pert_flow1.residual));

            pert_area2 = area;
            pert_flow2 = flow_data;
            pert_area2.at(Si) = area.at(Si) - pertA;

            getDomainResi(flo_opts, pert_area2, pert_flow2.W, &(pert_flow2.fluxes), &(pert_flow2.residual));

            for (int k = 0; k < 3; k++) {
                const int ki = Ri * 3 + k;
                const int ki_glob = (Ri-1) * 3 + k;
                dRdArea(ki_glob, Si) = (pert_flow1.residual.at(ki) - pert_flow2.residual.at(ki)) / (2 * pertA);
            }
        }
    }

    return dRdArea;
}
SparseMatrix<double> evaldRdW_FD(
    const std::vector<double> &area,
	const struct Flow_options &flo_opts,
	const struct Flow_data &flow_data)
{
	const int n_elem = flo_opts.n_elem;
    const int n_resi = n_elem*3;
	const int n_face = n_elem+1;
    const int n_flux = n_face*3;

    SparseMatrix<double> dRdW(n_resi, n_resi);
    int Ri, Wi;
    int irow_glob, icol_glob;
    struct Flow_data pert_flow1 = flow_data;
    struct Flow_data pert_flow2 = flow_data;
    std::vector<double> dRdW_block(9, 0);
    int ki, kip;
    const double dh = 1e-07;
    // DR/DW
    for (int Ri = 1; Ri < n_elem+1; Ri++) {
        for (int Wi = 1; Wi < n_elem+1; Wi++) {
            for (int istate_w = 0; istate_w < 3; istate_w++) {

                const int first_cell = 1, last_cell = n_elem;
                const double dx = 1.0/n_elem;

                const double pertW = flow_data.W.at(Wi * 3 + istate_w) * dh;

                pert_flow1.W = flow_data.W;
                pert_flow1.W.at(Wi * 3 + istate_w) = flow_data.W.at(Wi * 3 + istate_w) + pertW;

                for (int i = 1; i < n_elem+1; i++) {
                    const double u = pert_flow1.W[i*3+1] / pert_flow1.W[i*3+0];
                    const double c = get_c(flo_opts.gam, pert_flow1.W[i*3+0], pert_flow1.W[i*3+1], pert_flow1.W[i*3+2]);
                    pert_flow1.dt[i] = (flo_opts.CFL * dx) / fabs(u + c);
                }
                pert_flow1.dt[first_cell] = pert_flow1.dt[first_cell+1];
                pert_flow1.dt[last_cell] = pert_flow1.dt[last_cell-1];

                for (int i = 0; i < 1000; i++) {
                    inletBC(flo_opts, pert_flow1.dt[first_cell], dx, &pert_flow1);
                    outletBC(flo_opts, pert_flow1.dt[last_cell], dx, &pert_flow1);
                }
                getDomainResi(flo_opts, area, pert_flow1.W, &(pert_flow1.fluxes), &(pert_flow1.residual));

                pert_flow2.W = flow_data.W;
                pert_flow2.W.at(Wi * 3 + istate_w) = flow_data.W.at(Wi * 3 + istate_w) - pertW;
                for (int i = 1; i < n_elem+1; i++) {
                    const double u = pert_flow2.W[i*3+1] / pert_flow2.W[i*3+0];
                    const double c = get_c(flo_opts.gam, pert_flow2.W[i*3+0], pert_flow2.W[i*3+1], pert_flow2.W[i*3+2]);
                    pert_flow2.dt[i] = (flo_opts.CFL * dx) / fabs(u + c);
                }
                pert_flow2.dt[first_cell] = pert_flow2.dt[first_cell+1];
                pert_flow2.dt[last_cell] = pert_flow2.dt[last_cell-1];

                for (int i = 0; i < 1000; i++) {
                    inletBC(flo_opts, pert_flow2.dt[first_cell], dx, &pert_flow2);
                    outletBC(flo_opts, pert_flow2.dt[last_cell], dx, &pert_flow2);
                }
                getDomainResi(flo_opts, area, pert_flow2.W, &(pert_flow2.fluxes), &(pert_flow2.residual));

                for (int istate_resi = 0; istate_resi < 3; istate_resi++) {
                    ki = Ri * 3 + istate_resi;
                    dRdW_block.at(istate_resi * 3 + istate_w) = (pert_flow1.residual.at(ki) - pert_flow2.residual.at(ki)) / (2 * pertW);
                }

            } // END STATEI LOOP

            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    irow_glob = (Ri-1) * 3 + row;
                    icol_glob = (Wi-1) * 3 + col;
                    dRdW.insert(irow_glob, icol_glob) = dRdW_block.at(row * 3 + col);
                }
            }
        }
    }
    return dRdW;
}
