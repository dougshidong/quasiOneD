#include<Eigen/Core>
#include<Eigen/Sparse>
#include<math.h>
#include<iostream>
#include"globals.h"
#include"convert.h"
#include"flux.h"
#include"quasiOneD.h"
#include"BCDerivatives.h"

using namespace Eigen;

// Forward Declarations
void ScalarJac(
    std::vector <double> W,
    std::vector <double> &Ap_list,
    std::vector <double> &An_list);

void BCJac(
    std::vector <double> W,
    std::vector <double> dt,
    std::vector <double> dx,
    std::vector <double> &dBidWi,
    std::vector <double> &dBidWd,
    std::vector <double> &dBodWd,
    std::vector <double> &dBodWo);

std::vector <double> evaldQdW(
    std::vector <double> W,
    std::vector <double> S);

std::vector <double> evaldpdW(
    std::vector <double> W,
    std::vector <double> S);

//************************************************************************
// Calculates Jacobian
SparseMatrix<double> evaldRdW(
    std::vector <double> W,
    std::vector <double> dx,
    std::vector <double> dt,
    std::vector <double> S)
{
    SparseMatrix<double> dRdW(3 * nx, 3 * nx);
    dRdW.reserve(27 * (nx - 2) + 27 * 4);
    double rhoin = W[0];
    double uin = W[1] / rhoin;
    double ein = W[2];
    double pin = (gam - 1) * ( ein - rhoin * uin * uin / 2.0 );
    double cin = sqrt( gam * pin / rhoin );
    double Min = uin/cin;
    int Ri, Wi;
    int k, rowi, coli;
    double val;
    // Get Jacobians and Fluxes
    std::vector <double> Ap(nx * 3 * 3, 0), An(nx * 3 * 3, 0);
    if(FluxScheme == 0) ScalarJac(W, Ap, An);
    // Get Boundary Jacobians
    // DO NOT USE IMPLICIT FLOW SOLVERS WITH ADJOINT
    std::vector <double> dBidWi(9), dBidWd(9), dBodWd(9), dBodWo(9);
    if(StepScheme < 3)  BCJac(W, dt, dx, dBidWi, dBidWd, dBodWd, dBodWo);
    // Evaluate dpdW
    std::vector <double> dpdW(3 * nx, 0);
    dpdW = evaldpdW(W, S);
    // Input 4 lines where BC Jacobians occur
    // psi(1), psi(2), psi(n-1), psi(n)
    for(int row = 0; row < 3; row++)
    {
        for(int col = 0; col < 3; col++)
        {
            k = row * 3 + col;
            // d(inlet)/d(inlet)
            // R0, W0
            Ri = 0;
            Wi = 0;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;

            val = - dBidWi[k];
            dRdW.insert(rowi, coli) = val;

            // d(inlet)/d(domain)
            // R0, W1
            Ri = 0;
            Wi = 1;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;

            val = - dBidWd[k];
            dRdW.insert(rowi, coli) = val;

            // d(outlet)/d(outlet)
            // R = nx - 1, W = nx - 1
            Ri = nx - 1;
            Wi = nx - 1;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;

            val = - dBodWo[k];
            dRdW.insert(rowi, coli) = val;

            // d(outlet)/d(domain)
            // R = nx - 1, W = nx - 2
            Ri = nx - 1;
            Wi = nx - 2;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;

            val = - dBodWd[k];
            dRdW.insert(rowi, coli) = val;
        }
    }
    for(int Ri = 1; Ri < nx - 1; Ri++)
    {
        Wi = Ri - 1;
        if(Wi >= 0)
        {
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                k = row * 3 + col;
                rowi = Ri * 3 + row;
                coli = Wi * 3 + col;

                val = - Ap[Wi * 9 + k] * S[Ri];
                dRdW.insert(rowi, coli) = val;
            }
        }

        Wi = Ri;
        if(Wi >= 0 && Wi <= nx - 1)
        {
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                k = row * 3 + col;
                rowi = Ri * 3 + row;
                coli = Wi * 3 + col;

                val = Ap[Wi * 9 + k] * S[Ri + 1];
                val -= An[Wi * 9 + k] * S[Ri];
                if(row == 1)
                {
                    val -= dpdW[Wi * 3 + col] * (S[Ri + 1] - S[Ri]);
                }

                dRdW.insert(rowi, coli) = val;
            }
        }

        Wi = Ri + 1;
        if(Wi <= nx - 1)
        {
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                k = row * 3 + col;
                rowi = Ri * 3 + row;
                coli = Wi * 3 + col;

                val = An[Wi * 9 + k] * S[Ri + 1];
                dRdW.insert(rowi, coli) = val;
            }
        }
    }
    if(Min > 1.0)
    {
        // Supersonic Inlet, don't solve for psi(0)
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            // R1, W0
            Ri = 1;
            Wi = 0;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;

            dRdW.coeffRef(rowi, coli) = 0;
        }
    }
    return dRdW;
}

SparseMatrix<double> evaldRdW_FD(
    std::vector <double> W,
    std::vector <double> S)
{
    SparseMatrix<double> dRdW(3 * nx, 3 * nx);
    double rhoin = W[0];
    double uin = W[1] / rhoin;
    double ein = W[2];
    double pin = (gam - 1) * ( ein - rhoin * uin * uin / 2.0 );
    double cin = sqrt( gam * pin / rhoin );
    double Min = uin/cin;
    int Ri, Wi;
    int rowi, coli;
    std::vector <double> Wd(3 * nx, 0), Q(3 * nx, 0);
    std::vector <double> Flux(3 * (nx + 1), 0);
    std::vector <double> Resi1(3 * nx, 0), Resi2(3 * nx, 0);
    std::vector <double> dRdW_block(9, 0), dRdWp(9, 0);
    WtoQ(W, Q, S);
    getFlux(Flux, W);
    int ki, kip;
    double pert;

    // DR/DW
    for(int Ri = 0; Ri < nx; Ri++) // LOOP OVER R
    {
        for(int Wi = 0; Wi < nx; Wi++) // LOOP OVER W
        {
            double h = 1e-8;
            for(int statei = 0; statei < 3; statei++) // LOOP OVER STATEI
            {
                for(int i = 0; i < 3 * nx; i++)
                    Wd[i] = W[i];

                pert = W[Wi * 3 + statei] * h;
                Wd[Wi * 3 + statei] = W[Wi * 3 + statei] + pert;

                // RESI 1
                // Inlet
                if(Ri == 0) inletBC(Wd, Resi1, 1.0, 1.0);
                // Outlet
                else if (Ri == nx - 1) outletBC(Wd, Resi1, 1.0, 1.0);
                // Domain
                else
                {
                    WtoQ(Wd, Q, S);
                    getFlux(Flux, Wd);

                    for(int resii = 0; resii < 3; resii++)
                    {
                        ki = Ri * 3 + resii;
                        kip = (Ri + 1) * 3 + resii;
                        Resi1[ki] = Flux[kip] * S[Ri + 1] - Flux[ki] * S[Ri] - Q[ki];
                    }
                }

                for(int i = 0; i < 3 * nx; i++)
                    Wd[i] = W[i];

                Wd[Wi * 3 + statei] = W[Wi * 3 + statei] - pert;
                // RESI 2
                // Inlet
                if(Ri == 0) inletBC(Wd, Resi2, 1.0, 1.0);
                // Outlet
                else if (Ri == nx - 1) outletBC(Wd, Resi2, 1.0, 1.0);
                // Domain
                else
                {
                    WtoQ(Wd, Q, S);
                    getFlux(Flux, Wd);

                    for(int resii = 0; resii < 3; resii++)
                    {
                        ki = Ri * 3 + resii;
                        kip = (Ri + 1) * 3 + resii;
                        Resi2[ki] = Flux[kip] * S[Ri + 1] - Flux[ki] * S[Ri] - Q[ki];
                    }
                }

                for(int resii = 0; resii < 3; resii++)
                {
                    ki = Ri * 3 + resii;
                    dRdW_block[resii * 3 + statei] = (Resi1[ki] - Resi2[ki]) / (2 * pert);
                }

            } // END STATEI LOOP

            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                rowi = Ri * 3 + row;
                coli = Wi * 3 + col;
                dRdW.insert(rowi, coli) = dRdW_block[row * 3 + col];
            }
        }  // END LOOP OVER W
    } // END LOOP OVER R
    if(Min > 1.0)
    {
        // Supersonic Inlet, don't solve for psi(0)
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            // R0, W0
            Ri = 0;
            Wi = 0;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;
            dRdW.coeffRef(rowi, coli) = 0;
            if(row == col)
                dRdW.coeffRef(rowi, coli) = 1;

            // R1, W0
            Ri = 1;
            Wi = 0;
            rowi = Ri * 3 + row;
            coli = Wi * 3 + col;

            dRdW.coeffRef(rowi, coli) = 0;
        }
    }
    return dRdW;
}

// Steger-Warming Flux Splitting
void StegerJac(
    std::vector <double> W,
    std::vector <double> &Ap_list,
    std::vector <double> &An_list,
    std::vector <double> &Flux)
{
    double eps = 0.1;
    double gam = 1.4;
    double M[3][3] = {{0}},
           Minv[3][3] = {{0}},
           N[3][3] = {{0}},
           Ninv[3][3] = {{0}},
           lambdaP[3][3],
           lambdaN[3][3];
    double lambdaa[3];


    double Ap[3][3], An[3][3], tempP[3][3], tempN[3][3], prefix[3][3], suffix[3][3];

    std::vector <double> rho(nx), u(nx), p(nx), c(nx);
    std::vector <double> Ap_list1(nx * 3 * 3, 0), An_list1(nx * 3 * 3, 0);

    double beta = gam - 1;

    for(int i = 0; i < nx; i++)
    {
        rho[i] = W[i * 3 + 0];
        u[i] = W[i * 3 + 1] / rho[i];
        p[i] = (gam-1) * (W[i * 3 + 2] - rho[i] * pow(u[i], 2) / 2);
        c[i] = sqrt( gam * p[i] / rho[i] );
    }


    for(int i = 0; i < nx; i++)
    {
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            Ap[row][col] = 0;
            An[row][col] = 0;
            tempP[row][col] = 0;
            tempN[row][col] = 0;
            prefix[row][col] = 0;
            suffix[row][col] = 0;
            lambdaP[row][col] = 0;
            lambdaN[row][col] = 0;
        }

        M[0][0] = 1.0;
        M[1][0] = -u[i] / rho[i];
        M[2][0] = 0.5 * u[i] * u[i] * beta;
        M[1][1] = 1.0 / rho[i];
        M[2][1] = -u[i] * beta;
        M[2][2] = beta;
        Minv[0][0] = 1.0;
        Minv[1][0] = u[i];
        Minv[2][0] = 0.5 * u[i] * u[i];
        Minv[1][1] = rho[i];
        Minv[2][1] = u[i] * rho[i];
        Minv[2][2] = 1.0 / beta;
        N[0][0] = 1.0;
        N[1][1] = rho[i] * c[i];
        N[2][1] = -rho[i] * c[i];
        N[0][2] = -1.0 / (c[i] * c[i]);
        N[1][2] = 1.0;
        N[2][2] = 1.0;
        Ninv[0][0] = 1.0;
        Ninv[0][1] = 1.0 / (2.0 * c[i] * c[i]);
        Ninv[0][2] = 1.0 / (2.0 * c[i] * c[i]);
        Ninv[1][1] = 1.0 / (2.0 * rho[i] * c[i]);
        Ninv[1][2] = -1.0 / (2.0 * rho[i] * c[i]);
        Ninv[2][1] = 0.5;
        Ninv[2][2] = 0.5;
        lambdaa[0] = u[i];
        lambdaa[1] = u[i] + c[i];
        lambdaa[2] = u[i] - c[i];

        for(int k = 0; k < 3; k++)
            if(lambdaa[k] > 0)
                lambdaP[k][k] = (lambdaa[k] + sqrt(pow(lambdaa[k], 2) + pow(eps, 2))) / 2.0;
            else
                lambdaN[k][k] = (lambdaa[k] - sqrt(pow(lambdaa[k], 2) + pow(eps, 2))) / 2.0;

        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                prefix[row][col]+= Minv[row][k] * Ninv[k][col];
                suffix[row][col]+= N[row][k] * M[k][col];
            }
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                tempP[row][col] += prefix[row][k] * lambdaP[k][col];
                tempN[row][col] += prefix[row][k] * lambdaN[k][col];
            }
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                Ap[row][col]+= tempP[row][k] * suffix[k][col];
                An[row][col]+= tempN[row][k] * suffix[k][col];
            }
        // could remove above loop and just use aplist and anlist
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            int vec_pos = (i * 3 * 3) + (row * 3) + col;
            Ap_list1[vec_pos] = Ap[row][col];
            An_list1[vec_pos] = An[row][col];
        }

    }

    for(int i = 1; i < nx; i++)
    {
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            int Ap_pos = ((i - 1) * 3 * 3) + (row * 3) + col;
            int An_pos = (i * 3 * 3) + (row * 3) + col;
            Flux[i * 3 + row] += Ap_list1[Ap_pos] * W[(i - 1) * 3 + col]
                                 + An_list1[An_pos] * W[i * 3 + col];
        }
    }

}

void JacobianCenter(
    std::vector <double> &J,
    double u, double c)
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

std::vector <double> evaldlambdadW(
    std::vector <double> W,
    int i)
{
    std::vector <double> dlambdadWp(3), dlambdadW(3);
    double rho, u, e, p, c;
    rho = W[i * 3 + 0];
    u = W[i * 3 + 1] / rho;
    e = W[i * 3 + 2];
    p = (gam - 1) * ( e - rho * u * u / 2 );
    c = sqrt( gam * p / rho );

    // dlambda/dWp
    double dlambdadr, dlambdadu, dlambdadp;

    dlambdadr = -c / (4.0 * rho);
    dlambdadu = 0.5;
    dlambdadp = c / (4.0 * p);

    dlambdadWp[0] = dlambdadr;
    dlambdadWp[1] = dlambdadu;
    dlambdadWp[2] = dlambdadp;

    std::vector <double> dwpdw(9, 0);
    dWpdW(dwpdw, W, i);
    dlambdadW[0] = 0;
    dlambdadW[1] = 0;
    dlambdadW[2] = 0;
    for(int row = 0; row < 1; row++)
    for(int col = 0; col < 3; col++)
    for(int k = 0; k < 3; k++)
        dlambdadW[row * 3 + col] += dlambdadWp[row * 3 + k] * dwpdw[k * 3 + col];

    return dlambdadW;
}

void ScalarJac(
    std::vector <double> W,
    std::vector <double> &Ap_list,
    std::vector <double> &An_list)
{
    std::vector <double> rho(nx), u(nx), e(nx);
    std::vector <double> T(nx), p(nx), c(nx), Mach(nx);
    WtoP(W, rho, u, e, p, c, T);

    int vec_pos, k;
    double lamb;

    std::vector <double> J(9, 0);
    std::vector <double> dlambdaPdW(3, 0);
    std::vector <double> dlambdaNdW(3, 0);
    // A+
    for(int i = 0; i < nx - 1; i++)
    {
        // dF/dW
        JacobianCenter(J, u[i], c[i]);

        // lambda
        lamb = (u[i] + u[i + 1] + c[i] + c[i + 1]) / 2.0;
//      // dlambdaP/dW
        dlambdaPdW = evaldlambdadW(W, i);

        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            vec_pos = (i * 9) + (row * 3) + col; // NOT Transposed
            k = row * 3 + col;
            Ap_list[vec_pos] = J[k] / 2.0 - dlambdaPdW[col] * Scalareps
                               * (W[(i + 1) * 3 + row] - W[i * 3 + row]) / 2.0;
            if(row == col)
            {
                Ap_list[vec_pos] += Scalareps * lamb / 2.0;
            }
        }
    }

    // A-
    for(int i = 1; i < nx; i++)
    {
        // dF/dW
        JacobianCenter(J, u[i], c[i]);

        // lambda
        lamb = (u[i] + u[i - 1] + c[i] + c[i - 1]) / 2.0;
//      // dlambdaN/dW
        dlambdaNdW = evaldlambdadW(W, i);

        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            vec_pos = (i * 9) + (row * 3) + col; // NOT Transposed
            k = row * 3 + col;
            An_list[vec_pos] = J[k] / 2.0 - dlambdaNdW[col] * Scalareps
                               * (W[i * 3 + row] - W[(i - 1) * 3 + row]) / 2.0;
            if(row == col)
            {
                An_list[vec_pos] -= Scalareps * lamb / 2.0;
            }
        }
    }
}

std::vector <double> evaldQdW(
    std::vector <double> W,
    std::vector <double> S)
{
    double dpdw[3], rho, u, dS;
    std::vector <double> dQdW(3 * nx);
    for(int i = 0; i < nx; i++)
    {
        rho = W[i * 3 + 0];
        u = W[i * 3 + 1] / rho;

        dpdw[0] = (gam - 1) / 2.0 * u * u;
        dpdw[1] = - (gam - 1) * u;
        dpdw[2] = (gam - 1);

        dS = S[i + 1] - S[i];

        dQdW[i * 3 + 0] = dpdw[0] * dS;
        dQdW[i * 3 + 1] = dpdw[1] * dS;
        dQdW[i * 3 + 2] = dpdw[2] * dS;
    }
    return dQdW;
}

std::vector <double> evaldpdW(
    std::vector <double> W,
    std::vector <double> S)
{
    std::vector <double> dpdW(3 * nx);
    double rho, u;
    for(int i = 0; i < nx; i++)
    {
        rho = W[i * 3 + 0];
        u = W[i * 3 + 1] / rho;

        dpdW[i * 3 + 0] = (gam - 1) / 2.0 * u * u;
        dpdW[i * 3 + 1] = - (gam - 1) * u;
        dpdW[i * 3 + 2] = (gam - 1);
    }
    return dpdW;
}
// DERIVATIVES WRT GEOMETRY

MatrixXd evaldRdS(
    std::vector <double> Flux,
    std::vector <double> S,
    std::vector <double> W)
{
    MatrixXd dRdS(3 * nx, nx + 1);
    std::vector <double> Q(3 * nx, 0), p(nx);
    WtoQ(W, Q, S);
    getp(W, p);
    int Si, kR, kS;
    dRdS.setZero();
    for(int Ri = 1; Ri < nx - 1; Ri++)
    {
        for(int k = 0; k < 3; k++)
        {
            kR = Ri * 3 + k;

            Si = Ri;
            kS = Si * 3 + k;
            dRdS(kR, Si) = -Flux[kS];
            if(k == 1) dRdS(kR, Si) += p[Ri];

            Si = Ri + 1;
            kS = Si * 3 + k;
            dRdS(kR, Si) = Flux[kS];
            if(k == 1) dRdS(kR, Si) += -p[Ri];
        }
    }
    return dRdS;
}

MatrixXd evaldRdS_FD(
    std::vector <double> Flux,
    std::vector <double> S,
    std::vector <double> W)
{
    MatrixXd dRdS(3 * nx, nx + 1);
    dRdS.setZero();
    std::vector <double> Resi0(3 * nx, 0), Resi1(3 * nx, 0), Resi2(3 * nx, 0);
    std::vector <double> Sd(nx + 1, 0);
    std::vector <double> Q(3 * nx, 0);
    double h = 0.000000001;
    double pert;
    int ki, kip;
    for(int Ri = 1; Ri < nx - 1; Ri++)
    {
        for(int Si = 0; Si < nx + 1; Si++)
        {
            for(int m = 0; m < nx + 1; m++)
                Sd[m] = S[m];

            pert = S[Si] * h;
            Sd[Si] = S[Si] + pert;

            WtoQ(W, Q, Sd);

            for(int k = 0; k < 3; k++)
            {
                ki = Ri * 3 + k;
                kip = (Ri + 1) * 3 + k;
                Resi1[ki] = Flux[kip] * Sd[Ri + 1] - Flux[ki] * Sd[Ri] - Q[ki];
            }

            for(int m = 0; m < nx + 1; m++)
                Sd[m] = S[m];

            Sd[Si] = S[Si] - pert;

            WtoQ(W, Q, Sd);

            for(int k = 0; k < 3; k++)
            {
                ki = Ri * 3 + k;
                kip = (Ri + 1) * 3 + k;
                Resi2[ki] = Flux[kip] * Sd[Ri + 1] - Flux[ki] * Sd[Ri] - Q[ki];
                dRdS(ki, Si) = (Resi1[ki] - Resi2[ki]) / (2 * pert);
            }
        }
    }

    return dRdS;
}
