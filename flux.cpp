#include "flux.hpp"
#include "structures.hpp"
#include<vector>
#include<math.h>
#include<iostream>
#include"flux.hpp"
#include"convert.hpp"
#include<adolc/adolc.h>
#include <complex>

template<typename dreal>
void matrixMult(dreal A[][3], dreal B[][3], dreal result[][3]);

template<typename dreal>
void Flux_Scalar(
	const struct Flow_options &flow_options,
    const std::vector<dreal> &W,
    std::vector<dreal> *const fluxes);

// Get fluxes based on flux_scheme
template<typename dreal>
void getFlux(
	const struct Flow_options &flow_options,
	const std::vector<dreal> &W,
	std::vector<dreal> *const fluxes)
{
    if (flow_options.flux_scheme == 0)      // Scalar
        Flux_Scalar(flow_options, W, fluxes);
	else abort();
    //else if (flux_scheme == 1) // SW
    //    Flux_SW(fluxes, W);
    //else if (flux_scheme == 2) // MWS
    //    Flux_MSW(fluxes, W);
    //else if (flux_scheme == 3) // CMWS
    //    Flux_CMSW(fluxes, W);
    //else if (flux_scheme == 4) // Roe
    //    Flux_Roe(fluxes, W);
}
template void getFlux( const struct Flow_options &flow_options, const std::vector<double> &W, std::vector<double> *const fluxes);
template void getFlux( const struct Flow_options &flow_options, const std::vector<adouble> &W, std::vector<adouble> *const fluxes);
template void getFlux( const struct Flow_options &flow_options, const std::vector<std::complex<double>> &W, std::vector<std::complex<double>> *const fluxes);

template<typename dreal>
void Flux_Scalar(
	const struct Flow_options &flow_options,
    const std::vector<dreal> &W,
    std::vector<dreal> *const fluxes)
{
	const int n_face = flow_options.n_elem+1;
    assert(n_face == (*fluxes).size()/3);

    std::vector<dreal> F_m(3);
    std::vector<dreal> F_p(3);

	double gam = flow_options.gam;
	dreal w1 = W[0*3+0], w2 = W[0*3+1], w3 = W[0*3+2];
	dreal u_m = w2 / w1;
	dreal c_m = get_c(gam,w1,w2,w3);
	F_m[0] = w2;
	F_m[1] = w2 * w2 / w1 + (gam - 1.0) * ( w3 - w2 * w2 / (2.0 * w1) );
	F_m[2] = ( w3 + (gam - 1.0) * (w3 - w2 * w2 / (2.0 * w1)) ) * w2 / w1;

    for (int i_face = 0; i_face < n_face; i_face++) {
        const int i_cell = i_face+1;
		w1 = W[i_cell*3+0];
		w2 = W[i_cell*3+1];
		w3 = W[i_cell*3+2];
		dreal u_p = w2 / w1;
		dreal c_p = get_c(gam,w1,w2,w3);

        dreal avgu = ( u_m+u_p ) / 2.0;
        dreal avgc = ( c_m+c_p ) / 2.0;
        //lambda = std::max( std::max( fabs(avgu), fabs(avgu + avgc) ), fabs(avgu - avgc) );
        dreal lambda = avgu + avgc;

		F_p[0] = w2;
		F_p[1] = w2 * w2 / w1 + (gam - 1.0) * ( w3 - w2 * w2 / (2.0 * w1) );
		F_p[2] = ( w3 + (gam - 1.0) * (w3 - w2 * w2 / (2.0 * w1)) ) * w2 / w1;

        for (int i_state = 0; i_state < 3; i_state++) {
            const int kip = (i_face + 1) * 3 + i_state;
            const int ki  = (i_face - 0) * 3 + i_state;
            (*fluxes)[ki] = 0.5 * (F_m[i_state] + F_p[i_state]) - 0.5 * flow_options.scalar_d_eps * lambda * (W[kip] - W[ki]);
        }
		F_m = F_p;
		u_m = u_p;
		c_m = c_p;
    }
    return;
}

//void initializeFlux(int n_elem)
//{
//    if (flux_scheme == 0)
//    {
//        u.resize(n_elem);
//        c.resize(n_elem);
//        F.resize(3 * n_elem);
//    }
//  if (flux_scheme == 1 || flux_scheme == 2 || flux_scheme == 3)
//  {
//      Ap_list.resize(n_elem * 3 * 3);
//      An_list.resize(n_elem * 3 * 3);
//      Wavg.resize(3 * n_elem);
//  }
//  if (flux_scheme == 3)
//  {
//      FluxSW.resize(3 * (n_elem + 1));
//      FluxMSW.resize(3 * (n_elem + 1));
//      p.resize(n_elem);
//  }
//  if (flux_scheme == 4)
//  {
//      F.resize(3 * n_elem);
//  }
//}

//void evalLambda(dreal lamb[], dreal u, dreal c)
//{
//    lamb[0] = u;
//    lamb[1] = u + c;
//    lamb[2] = u - c;
//
//    return;
//}
//
//// fluxes Jacobian for SW/MSW/CMSW
//void Flux_Jacobian(
//    std::vector<dreal> &Ap_list,
//    std::vector<dreal> &An_list,
//    std::vector<dreal> const &W)
//{
//    dreal eps = 0.1;
//
//    dreal S[3][3] = {{0}},
//           Sinv[3][3] = {{0}},
//           C[3][3] = {{0}},
//           Cinv[3][3] = {{0}},
//           lambdaP[3][3],
//           lambdaN[3][3];
//    dreal lambdaa[3];
//
//
//    dreal Ap[3][3], An[3][3], tempP[3][3], tempN[3][3], prefix[3][3], suffix[3][3];
//
//    dreal beta = 0.4;//gam-1;
//    dreal rho, u, e, c;
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
//        rho = W[i * 3 + 0];
//        u = W[i * 3 + 1] / rho;
//        e = W[i * 3 + 2];
//        c = sqrt( gam / rho * (gam - 1.0) * ( e - rho * u * u / 2 ) );
//        S[0][0] = 1.0;
//        S[1][0] = -u / rho;
//        S[2][0] = 0.5 * u * u * beta;
//        S[1][1] = 1.0 / rho;
//        S[2][1] = -u * beta;
//        S[2][2] = beta;
//        Sinv[0][0] = 1.0;
//        Sinv[1][0] = u;
//        Sinv[2][0] = 0.5 * u * u;
//        Sinv[1][1] = rho;
//        Sinv[2][1] = u * rho;
//        Sinv[2][2] = 1.0 / beta;
//        C[0][0] = 1.0;
//        C[1][1] = rho * c;
//        C[2][1] = -rho * c;
//        C[0][2] = -1.0 / (c * c);
//        C[1][2] = 1.0;
//        C[2][2] = 1.0;
//        Cinv[0][0] = 1.0;
//        Cinv[0][1] = 1.0 / (2.0 * c * c);
//        Cinv[0][2] = 1.0 / (2.0 * c * c);
//        Cinv[1][1] = 1.0 / (2.0 * rho * c);
//        Cinv[1][2] = -1.0 / (2.0 * rho * c);
//        Cinv[2][1] = 0.5;
//        Cinv[2][2] = 0.5;
//        evalLambda(lambdaa, u, c);
//
//        for (int k = 0; k < 3; k++)
//            if (lambdaa[k] > 0)
//                lambdaP[k][k] = (lambdaa[k] + sqrt(pow(lambdaa[k], 2) + pow(eps, 2))) /2;
//            else
//                lambdaN[k][k] = (lambdaa[k] - sqrt(pow(lambdaa[k], 2) + pow(eps, 2))) / 2;
//
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//            for (int k = 0; k < 3; k++)
//            {
//                prefix[row][col] += Sinv[row][k] * Cinv[k][col];
//                suffix[row][col] += C[row][k] * S[k][col];
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
//                Ap[row][col] += tempP[row][k] * suffix[k][col];
//                An[row][col] += tempN[row][k] * suffix[k][col];
//            }
//        // could remove above loop and just use aplist and anlist
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//        {
//            int vec_pos = (i * 3 * 3) + (row * 3) + col;
//            Ap_list[vec_pos] = Ap[row][col];
//            An_list[vec_pos] = An[row][col];
//        }
//
//    }
//}
//
//// StegerWarming
//void Flux_SW(
//    std::vector<dreal> &fluxes,
//    std::vector<dreal> const &W)
//{
//
//    Flux_Jacobian(Ap_list, An_list, W);
//
//    for (int i = 1; i < n_elem; i++)
//    {
//        fluxes[i * 3 + 0] = 0;
//        fluxes[i * 3 + 1] = 0;
//        fluxes[i * 3 + 2] = 0;
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//        {
//            int Ap_pos = ((i - 1) * 3 * 3) + (row * 3) + col;
//            int An_pos = (i * 3 * 3) + (row * 3) + col;
//            fluxes[i * 3 + row] += Ap_list[Ap_pos] * W[(i - 1) * 3 + col]
//                               + An_list[An_pos] * W[i * 3 + col];
//        }
//    }
//}
//
//// Modified StegerWarming
//void Flux_MSW(
//    std::vector<dreal> &fluxes,
//    std::vector<dreal> const &W)
//{
//
//    for (int i = 0; i < n_elem - 1; i++)
//    {
//        for (int k = 0; k < 3; k++)
//        {
//            Wavg[i * 3 + k] = (W[i * 3 + k] + W[(i+1) * 3 + k]) / 2.0;
//        }
//    }
//    Wavg[(n_elem-1)*3 + 0] = 1;
//    Wavg[(n_elem-1)*3 + 1] = 1;
//    Wavg[(n_elem-1)*3 + 2] = 1;
//
//    Flux_Jacobian(Ap_list, An_list, Wavg);
//
//    for (int i = 1; i < n_elem; i++)
//    {
//        fluxes[i * 3 + 0] = 0;
//        fluxes[i * 3 + 1] = 0;
//        fluxes[i * 3 + 2] = 0;
//        for (int row = 0; row < 3; row++)
//        for (int col = 0; col < 3; col++)
//        {
//            int Ahalf_pos = ((i-1) * 3 * 3) + (row * 3) + col;
//            fluxes[i * 3 + row] += Ap_list[Ahalf_pos] * W[(i - 1) * 3 + col]
//                               + An_list[Ahalf_pos] * W[i * 3 + col];
//        }
//    }
//}
//
//// Corrected Modified StegerWarming
//void Flux_CMSW(
//    std::vector<dreal> &fluxes,
//    std::vector<dreal> const &W)
//{
//    Flux_SW(FluxSW, W);
//    Flux_MSW(FluxMSW, W);
//
//    for (int i = 0; i < n_elem; i++)
//    {
//        p[i] = (gam - 1.0) * ( W[i * 3 + 2] - (pow(W[i * 3 + 1], 2.0) / W[i * 3 + 0]) / 2.0 );
//    }
//
//    for (int i = 0; i < n_elem - 1; i++)
//    {
//        for (int k = 0; k < 3; k++)
//        {
//            Wavg[i * 3 + k]= 1.0 / (1.0 + pow( (p[i+1] - p[i]) / std::min(p[i+1], p[i]), 2 ));
//        }
//    }
//
//    for (int i = 1; i < n_elem; i++)
//    {
//        for (int k = 0; k < 3; k++)
//        {
//            fluxes[i * 3 + k] = 0;
//            fluxes[i * 3 + k] = Wavg[(i - 1) * 3 + k] * FluxMSW[i * 3 + k]
//                     + (1.0 - Wavg[(i - 1) * 3 + k]) * FluxSW[i * 3 + k];
//        }
//    }
//}
//
//
//
//void Flux_Roe(
//    std::vector<dreal> &fluxes,
//    std::vector<dreal> const &W)
//{
//    // Get Convective Variables
//    WtoF_all(W, F);
//
//    dreal beta = gam - 1.0;
//    dreal rH, uH, hH, cH;
//    dreal r1, u1, e1, p1, c1;
//    dreal r2, u2, e2, p2, c2;
//    dreal sr1, sr2;
//    dreal temp;
//    dreal epsilon[3];
//    dreal lamb[3], lambp1[3], lambH[3];
//    dreal lambdaP[3][3], lambdaN[3][3];
//    dreal S[3][3], Sinv[3][3], C[3][3], Cinv[3][3];
//    dreal Ap[3][3], An[3][3];
//    for (int row = 0; row < 3; row++)
//    {
//        for (int col = 0; col < 3; col++)
//        {
//            S[row][col]=0;
//            C[row][col]=0;
//            Cinv[row][col]=0;
//            Sinv[row][col]=0;
//            Ap[row][col]=0;
//            An[row][col]=0;
//            lambdaP[row][col] = 0;
//            lambdaN[row][col] = 0;
//        }
//    }
//    int i1, i2;
//    for (int i = 0; i < n_elem - 1; i++)
//    {
//        for (int row = 0; row < 3; row++)
//        {
//            for (int col = 0; col < 3; col++)
//            {
//                S[row][col]=0;
//                C[row][col]=0;
//                Cinv[row][col]=0;
//                Sinv[row][col]=0;
//            }
//        }
//        for (int k = 0; k < 3; k++)
//        {
//            lambdaP[k][k] = 0;
//            lambdaN[k][k] = 0;
//        }
//
//        i1 = i * 3;
//        i2 = (i + 1) * 3;
//
//        r1 = W[i1 + 0];
//        u1 = W[i1 + 1] / r1;
//        e1 = W[i1 + 2];
//        p1 = (gam - 1.0) * ( e1 - r1 * u1 * u1 / 2.0 );
//        c1 = sqrt( p1 * gam / r1 );
//
//        r2 = W[i2 + 0];
//        u2 = W[i2 + 1] / r2;
//        e2 = W[i2 + 2];
//        p2 = (gam - 1.0) * ( e2 - r2 * u2 * u2 / 2.0 );
//        c2 = sqrt( p2 * gam / r2 );
//
//        sr1 = sqrt(r1);
//        sr2 = sqrt(r2);
//
//        rH = sr1 * sr2;
//        uH = (sr1 * u1 + sr2 * u2) / (sr1 + sr2);
//        hH = (sr1 * (e1 + p1) / r1 + sr2 * (e2 + p2) / r2) / (sr1 + sr2);
//        cH = sqrt((gam - 1.0) * (hH - uH * uH / 2.0));
//
//        S[0][0] = 1.0;
//        S[1][0] = -uH / rH;
//        S[2][0] = 0.5 * uH * uH * beta;
//        S[1][1] = 1.0 / rH;
//        S[2][1] = -uH * beta;
//        S[2][2] = beta;
//        Sinv[0][0] = 1.0;
//        Sinv[1][0] = uH;
//        Sinv[2][0] = 0.5 * uH * uH;
//        Sinv[1][1] = rH;
//        Sinv[2][1] = uH * rH;
//        Sinv[2][2] = 1.0 / beta;
//        C[0][0] = 1.0;
//        C[1][1] = rH * cH;
//        C[2][1] = -rH * cH;
//        C[0][2] = -1.0 / (cH * cH);
//        C[1][2] = 1.0;
//        C[2][2] = 1.0;
//        Cinv[0][0] = 1.0;
//        Cinv[0][1] = 1.0 / (2.0 * cH * cH);
//        Cinv[0][2] = 1.0 / (2.0 * cH * cH);
//        Cinv[1][1] = 1.0 / (2.0 * rH * cH);
//        Cinv[1][2] = -1.0 / (2.0 * rH * cH);
//        Cinv[2][1] = 0.5;
//        Cinv[2][2] = 0.5;
//
//        evalLambda(lamb, u1, c1);
//        evalLambda(lambp1, u2, c2);
//        evalLambda(lambH, uH, cH);
//        for (int k = 0; k < 3; k++)
//        {
//            epsilon[k] = std::max(std::max(
//                            0.0,
//                            lambH[k] - lamb[k] ),
//                            lambp1[k] - lambH[k] );
//            if (fabs(lambH[k]) <= epsilon[k])
//            {
//                lambH[k] = 0.5 * (lambH[k] * lambH[k] / epsilon[k] + epsilon[k]);
//            }
//
//            if (lambH[k] > 0)
//            {
//                lambdaP[k][k] = lambH[k];
//            }
//            else
//            {
//                lambdaN[k][k] = lambH[k];
//            }
//        }
//
////      matrixMult(Sinv, Cinv, Ap);
////      matrixMult(Ap, lambdaP, Ap);
////      matrixMult(Ap, C, Ap);
////      matrixMult(Ap, S, Ap);
////
////      matrixMult(Sinv, Cinv, An);
////      matrixMult(An, lambdaN, An);
////      matrixMult(An, C, An);
////      matrixMult(An, S, An);
//
//        matrixMult(Sinv, Cinv, Cinv);  // Prefix
//        matrixMult(C, S, C);           // Suffix
//
//        matrixMult(Cinv, lambdaP, Ap);
//        matrixMult(Ap, C, Ap);
//
//        matrixMult(Cinv, lambdaN, An);
//        matrixMult(An, C, An);
//
//        for (int row = 0; row < 3; row++)
//        {
//            temp = 0;
//            for (int col = 0; col < 3; col++)
//            {
//                temp += (Ap[row][col] - An[row][col]) * (W[i2 + col] - W[i1 + col]);
//            }
//            fluxes[i2 + row] =
//                0.5 * (F[i1 + row] + F[i2 + row]) - 0.5 * temp;
//        }
//
//    }
//}
//
//void matrixMult(dreal A[][3], dreal B[][3], dreal result[][3])
//{
//    dreal temp[3][3];
//    for (int row=0;row<3;row++)
//        for (int col=0;col<3;col++)
//            temp[row][col]=0;
//
//    for (int row=0;row<3;row++)
//        for (int col=0;col<3;col++)
//        {
//            for (int k=0;k<3;k++)
//                temp[row][col]+=A[row][k]*B[k][col];
//        }
//    for (int row=0;row<3;row++)
//        for (int col=0;col<3;col++)
//            result[row][col]=temp[row][col];
//}
//
