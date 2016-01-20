// Calculates the discrete costate fluxes
#include<iostream>
#include<math.h>
#include"quasiOneD.h"
#include<vector>
#include<Eigen/Dense>
#include<Eigen/LU>
#include <stdio.h>
#include <iomanip>
#include "globals.h"
#include "flovar.h"
#include "convert.h"

double adjConv = 1e-7;
int adjMaxIt = 960000;
int adjPrintIt = 10000; 
int adjPrintConv = 1; // 0 to hide real-time adjConvergence

using namespace Eigen;

FullPivLU<MatrixXd> eigenFullPivLU(std::vector <double> matA, int m, int n)
{
    MatrixXd Ae(m, n);
    for(int i = 0; i < m; i++)
    for(int j = 0; j < n; j++)
        Ae(i, j) = matA[i * n + j];
    FullPivLU<MatrixXd> lu(Ae);
    return lu;
}

void StegerJac(std::vector <double> W,
               std::vector <double> &Ap_list,
               std::vector <double> &An_list);

void BCJac(std::vector <double> W,
           std::vector <double> &dBidWi,
           std::vector <double> &dBidWd,
           std::vector <double> &dBodWd,
           std::vector <double> &dBodWo);

void adjFlux(std::vector <double> S,
             std::vector <double> dx,
             std::vector <double> Ap_list,
             std::vector <double> An_list,
             std::vector <double> psi,
             std::vector <double> &psiF);

void adjointEuler(std::vector <double> S,
                  std::vector <double> dt,
                  std::vector <double> psiF,
                  std::vector <double> &Resi,
                  std::vector <double> &psi);

void adjointBC(std::vector <double> &psi, 
               std::vector <double> dBidWi,
               std::vector <double> dBodWo,
               std::vector <double> A1p,
               std::vector <double> Anp,
               std::vector <double> S);

void evaldIcdW(std::vector <double> &dIcdW,
               std::vector <double> W,
               std::vector <double> S);

void adjointSource(std::vector <double> &adjQ,
                   std::vector <double> W,
                   std::vector <double> S);


void adjoint(std::vector <double> x, 
                             std::vector <double> dx, 
                             std::vector <double> S,
                             std::vector <double> W,
                             std::vector <double> &psi)
{
    std::vector <double> dIcdW(3 * nx, 0);
    std::vector <double> adjQ(3 * nx, 0);
    std::vector <double> dBidWi(3 * 3, 0);
    std::vector <double> dBidWd(3 * 3, 0);
    std::vector <double> dBodWd(3 * 3, 0);
    std::vector <double> dBodWo(3 * 3, 0);
    std::vector <double> Resi(3 * nx, 0);
    std::vector <double> psiF(3 * (nx + 1), 0);
    std::vector <double> dt(nx);

    std::vector <int> itV(adjMaxIt/adjPrintIt);
    std::vector <double> normV(adjMaxIt/adjPrintIt);

    std::vector <double> Ap_list(nx * 3 * 3, 0), An_list(nx * 3 * 3, 0);

    double normR = 1.0;
    int iterations = 0;

    // Initialize Costate Vector
    for(int i = 0; i < 3 * nx; i++)
    {
        psi[i] = 1.0;
    }

    MatrixXd matA(3, 3);
    MatrixXd matB(3, 3);
    matA << 1, 2, 3, 4,5,6,7,8,9;
    std::vector <double> aa(3 * 3, 0);
    aa[0] = 1;
    aa[1] = 2;
    aa[2] = 3;
    aa[3] = 0;
    aa[4] = 1;
    aa[5] = 4;
    aa[6] = 5;
    aa[7] = 6;
    aa[8] = 0;
    int m = 3, n = 3;
    FullPivLU<MatrixXd> LU;
    LU = eigenFullPivLU(aa, m, n);
    std::cout << LU.matrixLU() <<std::endl;

    // Evaluate dIcdW
    evaldIcdW(dIcdW, W, S);

    // Adjoint Source Term
    adjointSource(adjQ, W, S);

    // Get Flux Jacobians
    StegerJac(W, Ap_list, An_list);
    // Boundary Flux Jacobians
    BCJac(W, dBidWi, dBidWd, dBodWd, dBodWo);

    std::vector <double> Anp(3 * 3, 0);
    std::vector <double> A1p(3 * 3, 0);
    adjointBC(psi, dBidWi, dBodWo, A1p, Anp, S);

    while(normR > adjConv && iterations < adjMaxIt)
    {
        iterations++;

        if(iterations % adjPrintIt == 0) 
        {
            if(adjPrintConv == 1)
            {
                std::cout<<"Iteration "<<iterations
                         <<"   NormR "<<std::setprecision(15)<<normR<<std::endl;
            }
            itV[iterations / adjPrintIt - 1] = iterations;
            normV[iterations / adjPrintIt - 1] = normR;
        }

        // CALCULATE TIME STEP
        for(int i = 0; i < nx; i++)
            dt[i] = 0.0001;


        adjFlux(S, dx, Ap_list, An_list, psi, psiF);
        adjointEuler(S, dt, psiF, Resi, psi);

        // Calculating the norm of the first costate residual
        normR = 0;
        for(int i = 0; i < nx; i++)
            normR = normR + Resi[2 * nx + i] * Resi[2 * nx + i];
        normR = sqrt(normR);
    }

    std::cout<<"Adjoint Iterations = "<<iterations<<"   Density Residual = "<<normR<<std::endl;
    
    FILE *Results;
    Results = fopen("Adjoint.dat", "w");
    fprintf(Results, "%d\n", nx);
    for(int k = 0; k < 3; k++)
    for(int i = 0; i < nx; i++)
        fprintf(Results, "%.15f\n", psi[k * nx + i]);

    fclose(Results);
}

void adjointEuler(std::vector <double> S,
                    std::vector <double> dt,
                    std::vector <double> psiF,
                    std::vector <double> &Resi,
                    std::vector <double> &psi)
{
    int ki;
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
            Resi[ki] = psiF[ki+1] * S[i+1]
                       - psiF[ki] * S[i];
        }
        Resi[k * nx + 0] = 0;
        Resi[k * nx + (nx - 1)]= 0;
    }
    for(int k = 0; k < 3; k++)
    for(int i = 1; i < nx - 1; i++)
    {
        ki = k * nx + i;
        psi[ki] = psi[ki] - (dt[i]) * Resi[ki];
    }
    return;
}

void adjFlux(
                   std::vector <double> S,
                   std::vector <double> dx,
                   std::vector <double> Ap_list,
                   std::vector <double> An_list,
                   std::vector <double> psi,
                   std::vector <double> &psiF)
{
    // Calculate psiF
    for(int i = 1; i < nx; i++)
    {
        psiF[0 * nx + i] = 0;
        psiF[1 * nx + i] = 0;
        psiF[2 * nx + i] = 0;
        for(int row = 0; row < 3; row++)
        {
            for(int col = 0; col < 3; col++)
            {
                int Ap_pos = ((i-1) * 3 * 3)+(row * 3)+col;
                int An_pos = (i * 3 * 3)+(row * 3)+col;
                psiF[row * nx + i] = psiF[row * nx + i]
                                + Ap_list[Ap_pos]
                                 * ( psi[col * nx + (i+1)] - psi[col * nx + i] )
                                + An_list[An_pos]
                                * ( psi[col * nx + i] - psi[col * nx + (i-1)] );
            }
        }
    }
}

// Calculates Jacobian
// Steger-Warming Flux Splitting
void StegerJac(std::vector <double> W,
               std::vector <double> &Ap_list,
               std::vector <double> &An_list)
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



    double beta = 0.4;//gam-1;

    for(int i = 0; i < nx; i++)
    {
        rho[i] = W[0 * nx + i];
        u[i] = W[1 * nx + i]/rho[i];
        p[i] = (gam-1) * (W[2 * nx + i]-rho[i] * pow(u[i], 2)/2);
        c[i] = sqrt((gam * p[i])/rho[i]);
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
    
        M[0][0] = 1;
        M[1][0] = -u[i]/rho[i];
        M[2][0] = 0.5 * u[i] * u[i] * beta;
        M[1][1] = 1/rho[i];
        M[2][1] = -u[i] * beta;
        M[2][2] = beta;
        Minv[0][0] = 1;
        Minv[1][0] = u[i];
        Minv[2][0] = 0.5 * u[i] * u[i];
        Minv[1][1] = rho[i];
        Minv[2][1] = u[i] * rho[i];
        Minv[2][2] = 1/beta;
        N[0][0] = 1;
        N[1][1] = rho[i] * c[i];
        N[2][1] = -rho[i] * c[i];
        N[0][2] = -1/(c[i] * c[i]);
        N[1][2] = 1;
        N[2][2] = 1;
        Ninv[0][0] = 1;
        Ninv[0][1] = 1/(2 * c[i] * c[i]);
        Ninv[0][2] = 1/(2 * c[i] * c[i]);
        Ninv[1][1] = 1/(2 * rho[i] * c[i]);
        Ninv[1][2] = -1/(2 * rho[i] * c[i]);
        Ninv[2][1] = 0.5;
        Ninv[2][2] = 0.5;
        lambdaa[0] = u[i];
        lambdaa[1] = u[i]+c[i];
        lambdaa[2] = u[i]-c[i];
        
        for(int k = 0; k < 3; k++)
            if(lambdaa[k]>0)
                lambdaP[k][k] = (lambdaa[k]+
                    sqrt(pow(lambdaa[k], 2)+pow(eps, 2)))/2;
            else
                lambdaN[k][k] = (lambdaa[k]-
                    sqrt(pow(lambdaa[k], 2)+pow(eps, 2)))/2;

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
            Ap_list[vec_pos] = Ap[row][col];
            An_list[vec_pos] = An[row][col];
        }

    }


}

// Set the adjoint BC for target pressure
void adjointBC(std::vector <double> &psi, 
               std::vector <double> dBidWi,
               std::vector <double> dBodWo,
               std::vector <double> A1p,
               std::vector <double> Anp,
               std::vector <double> S)
{
  //
  //psi[0 * nx] = 0;
  //psi[1 * nx - 1] = 0;
  //psi[1 * nx] = - (p0-pTarget[0]) * dx[0]/(S[1]-S[0]);
  //psi[2 * nx - 1] = - (pn-pTarget[nx - 1]) * dx[nx - 1]/(S[nx]-S[nx - 1]);
  //psi[2 * nx] = 0;
  //psi[3 * nx - 1] = 0;
}

void evaldIcdW(std::vector <double> &dIcdW,
               std::vector <double> W,
               std::vector <double> S)
{
    std::vector <double> ptarget(nx, 0);
    double dpdw[3], rho, u, p, dS;
    ioTargetPressure(-1, ptarget);
    for(int i = 0; i < nx; i++)
    {
        rho = W[0 * nx + i];
        u = W[1 * nx + i] / rho;
        p = (gam - 1) * ( W[2 * nx + i] - rho * u * u / 2 );

        dpdw[0] = (gam - 1) / 2 * u * u;
        dpdw[1] = - (gam - 1) * u;
        dpdw[2] = (gam - 1);

        dS = S[i + 1] - S[i];

        dIcdW[0 * nx + i] = (p - ptarget[i]) * dpdw[0] * dS;
        dIcdW[1 * nx + i] = (p - ptarget[i]) * dpdw[1] * dS;
        dIcdW[2 * nx + i] = (p - ptarget[i]) * dpdw[2] * dS;
    }
}
void adjointSource(std::vector <double> &adjQ,
                   std::vector <double> W,
                   std::vector <double> S)
{
    double dpdw[3], rho, u, dS;
    for(int i = 0; i < nx; i++)
    {
        rho = W[0 * nx + i];
        u = W[1 * nx + i] / rho;

        dpdw[0] = (gam - 1) / 2 * u * u;
        dpdw[1] = - (gam - 1) * u;
        dpdw[2] = (gam - 1);

        dS = S[i + 1] - S[i];

        adjQ[0 * nx + i] = dpdw[0] * dS;
        adjQ[1 * nx + i] = dpdw[1] * dS;
        adjQ[2 * nx + i] = dpdw[2] * dS;
    }
}

void BCJac(std::vector <double> W,
           std::vector <double> &dBidWi,
           std::vector <double> &dBidWd,
           std::vector <double> &dBodWd,
           std::vector <double> &dBodWo)
{
    std::vector <double> rho(nx), u(nx), e(nx), p(nx), c(nx), T(nx);
    std::vector <double> dbdwp(9,0), dwpdw(9);

    for(int i = 1; i < 9; i++)
    {
        dBidWi[i] = 0;
        dBidWd[i] = 0;
        dBodWd[i] = 0;
        dBodWo[i] = 0;
    }

    WtoP(W, rho, u, e, p, c, T);

    // ************************
    // OUTLET JACOBIANS
    // ************************

    double i1, i2;
    double r1, r2, p1, p2, u1, u2, c1, c2, t1;
    i1 = nx - 1;
    i2 = nx - 2;
    r1 = rho[i1];
    r2 = rho[i2];
    p1 = p[i1];
    p2 = p[i2];
    u1 = u[i1];
    u2 = u[i2];
    c1 = c[i1];
    c2 = c[i2];
    t1 = T[i1];

    // Shorthand
    double gamr, fu, drho, dp, du, cr, uu;
    drho = r1 - r2;
    dp = p1 - p2;
    du = u1 - u2;
    cr = r1 * c1;
    uu = u1 * u1;

    // Speed of Sound
    double dc1dr1, dc2dr2, dc1dp1, dc2dp2;
    dc1dr1 = - p1 * gam / (2 * cr * r1);
    dc2dr2 = - p2 * gam / (2 * c2 * r2 * r2);
    dc1dp1 = gam / (2 * cr);
    dc2dp2 = gam / (2 * c2 * r2);

    double eig1, eig2, eig3;
    double deig1du1, deig1du2;
    double deig2dr1, deig2du1, deig2dp1, deig2dr2, deig2du2, deig2dp2;
    double deig3dr1, deig3du1, deig3dp1, deig3dr2, deig3du2, deig3dp2;
    // Eigenvalue
    eig1 = (u1 + u2) / 2;
    eig2 = eig1 + (c1 + c2);
    eig3 = eig1 - (c1 + c2);

    deig1du1 = 1.0 / 2.0;
    deig1du2 = 1.0 / 2.0;

    deig2dr1 = dc1dr1 / 2;
    deig2du1 = deig1du1;
    deig2dp1 = dc1dp1 / 2;
    deig2dr2 = dc2dr2 / 2;
    deig2du2 = deig1du2;
    deig2dp2 = dc2dp2 / 2;

    deig3dr1 = - dc1dr1 / 2;
    deig3du1 = deig1du1;
    deig3dp1 = - dc1dp1 / 2;
    deig3dr2 = - dc2dr2 / 2;
    deig3du2 = deig1du2;
    deig3dp2 = - dc2dp2 / 2;

    // Riemann invariants
    double R1, R2, R3;
    double dR1dr1, dR1du1, dR1dp1, dR1dr2, dR1du2, dR1dp2;
    double dR2dr1, dR2du1, dR2dp1, dR2dr2, dR2du2, dR2dp2;
    double dR3dr1, dR3du1, dR3dp1, dR3dr2, dR3du2, dR3dp2;
    R1 = - eig1 * (drho - dp / (c1 * c1));
    R2 = - eig2 * (dp + cr * du);
    R3 = - eig3 * (dp - cr * du);

    dR1dr1 = - eig1 * (1 + 2 * dp * dc1dr1 / pow(c1, 3) );
    dR1du1 = deig1du1 * (dp - c1 * c1 * drho) / (c1 * c1);
    dR1dp1 = eig1 * (c1 - 2 * dp * dc1dp1) / pow(c1, 3);
    dR1dr2 = eig1;
    dR1du2 = deig1du2 * (dp - c1 * c1 * drho) / (c1 * c1);
    dR1dp2 = - eig1 / (c1 * c1);

    dR2dr1 = - du * eig2 * (c1 + r1 * dc1dr1) - (dp + cr * du) * deig2dr1;
    dR2du1 = - cr * eig2 - (dp + cr * du) * deig2du1;
    dR2dp1 = - eig2 * (1 + du * r1 * dc1dp1) - (dp + cr * du) * deig2dp1;
    dR2dr2 = - (dp + cr * du) * deig2dr2;
    dR2du2 = - cr * eig2 - (dp + cr * du) * deig2du2;
    dR2dp2 = eig2 - (dp + cr * du) * deig2dp2;

    dR3dr1 = c1 * du * eig3 + du * eig3 * r1 * dc1dr1 + (-dp + cr * du) * deig3dr1; 
    dR3du1 = - cr * eig3 + (-dp + cr * du) * deig3du1;
    dR3dp1 = - eig3 + du * eig3 * r1 * dc1dp1 + (-dp + cr * du) * deig3dp1;
    dR3dr2 = (-dp + cr * du) * deig3dr2;
    dR3du2 = - cr * eig3 + (-dp + cr * du) * deig3du2;
    dR3dp2 = eig3 + (-dp + cr * du) * deig3dp2;

    // dp1/dt
    double dp1dt;
    double dp1dtdr1, dp1dtdu1, dp1dtdp1;
    double dp1dtdr2, dp1dtdu2, dp1dtdp2;
    if(u1 < c1)
    {
        dp1dt = 0;
        dp1dtdr1 = 0;
        dp1dtdu1 = 0;
        dp1dtdp1 = 0;
        dp1dtdr2 = 0;
        dp1dtdu2 = 0;
        dp1dtdp2 = 0;
    }
    else
    {
        dp1dt = (R2 + R3) / 2;
        dp1dtdr1 = (dR2dr1 + dR3dr1) / 2; 
        dp1dtdu1 = (dR2du1 + dR3du1) / 2; 
        dp1dtdp1 = (dR2dp1 + dR3dp1) / 2; 
        dp1dtdr2 = (dR2dr2 + dR3dr2) / 2; 
        dp1dtdu2 = (dR2du2 + dR3du2) / 2; 
        dp1dtdp2 = (dR2dp2 + dR3dp2) / 2; 
    }

    // drho1/dt
    double dr1dt;
    double dr1dtdr1, dr1dtdu1, dr1dtdp1;
    double dr1dtdr2, dr1dtdu2, dr1dtdp2;
    dr1dt = R1 + dp1dt / (c1 * c1);

    dr1dtdr1 = dR1dr1 + (-2 * (R2 + R3) * dc1dr1 + c1 * (dR2dr1 + dR3dr1)) 
               / ( 2 * pow(c1, 3) );
    dr1dtdu1 = dR1du1 + (dR2du1 + dR3du1) / ( 2 * pow(c1, 2) );
    dr1dtdp1 = dR1dp1 + (-2 * (R2 + R3) * dc1dp1 + c1 * (dR2dp1 + dR3dp1)) 
               / ( 2 * pow(c1, 3) );
    dr1dtdr2 = dR1dr2 + (dR2dr2 + dR3dr2) / ( 2 * pow(c1, 2) );
    dr1dtdu2 = dR1du2 + (dR2du2 + dR3du2) / ( 2 * pow(c1, 2) );
    dr1dtdp2 = dR1dp2 + (dR2dp2 + dR3dp2) / ( 2 * pow(c1, 2) );

    // du1/dt
    double du1dt;
    double du1dtdr1, du1dtdu1, du1dtdp1;
    double du1dtdr2, du1dtdu2, du1dtdp2;
    du1dt = (R2 - dp1dt) / (cr);

    du1dtdr1 = ( (dp1dt - R2) * r1 * dc1dr1
               + c1 * (dp1dt - R2 - r1 * dp1dtdr1 + r1 * dR2dr1) )
               / (cr * cr);
    du1dtdu1 = (dR2du1 - dp1dtdu1) / cr;
    du1dtdp1 = ( (dp1dt - R2) * dc1dp1 + c1 * (dR2dp1 - dp1dtdp1) ) / (cr * c1);
    du1dtdr2 = (dR2dr2 - dp1dtdr2) / cr;
    du1dtdu2 = (dR2du2 - dp1dtdu2) / cr;
    du1dtdp2 = (dR2dp2 - dp1dtdp2) / cr;

    // d(ru)1/dt
    double dru1dtdr1, dru1dtdu1, dru1dtdp1;
    double dru1dtdr2, dru1dtdu2, dru1dtdp2;
    dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1;
    dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1;
    dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1;
    dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2;
    dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2;
    dru1dtdp2 = u1 * dr1dtdp2 + r1 * du1dtdp2;

    // de1/dt
    double de1dtdr1, de1dtdu1, de1dtdp1;
    double de1dtdr2, de1dtdu2, de1dtdp2;
    // de1dt = dp1dt * Cv / R + u1 + cr * du1dt + uu * dr1dt / 2;
    de1dtdr1 = du1dt * u1 + dp1dtdr1 * Cv / R + uu * dr1dtdr1 / 2 + cr * du1dtdr1;
    de1dtdu1 = dr1dt * u1 + du1dt * r1 + dp1dtdu1 * Cv / R
              + uu * dr1dtdu1 / 2 + cr * du1dtdu1;
    de1dtdp1 = dp1dtdp1 * Cv / R + uu * dr1dtdp1 / 2 + cr * du1dtdp1;
    de1dtdr2 = dp1dtdr2 * Cv / R + u1 + cr * du1dtdr2 + uu * dr1dtdr2 / 2;
    de1dtdu2 = dp1dtdu2 * Cv / R + u1 + cr * du1dtdu2 + uu * dr1dtdu2 / 2;
    de1dtdp2 = dp1dtdp2 * Cv / R + u1 + cr * du1dtdp2 + uu * dr1dtdp2 / 2;

    dbdwp[0] = dr1dtdr1;
    dbdwp[1] = dr1dtdu1;
    dbdwp[2] = dr1dtdp1;
    dbdwp[3] = dru1dtdr1;
    dbdwp[4] = dru1dtdu1;
    dbdwp[5] = dru1dtdp1;
    dbdwp[6] = de1dtdr1;
    dbdwp[7] = de1dtdu1;
    dbdwp[8] = de1dtdp1;

    // Get Transformation Matrix
    dWpdW(dwpdw, W, nx - 1);
    for(int row = 0; row < 3; row++)
    for(int col = 0; col < 3; col++)
        for(int k = 0; k < 3; k++)
        {
            dBodWo[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];
        }
    
    dbdwp[0] = dr1dtdr2;
    dbdwp[1] = dr1dtdu2;
    dbdwp[2] = dr1dtdp2;
    dbdwp[3] = dru1dtdr2;
    dbdwp[4] = dru1dtdu2;
    dbdwp[5] = dru1dtdp2;
    dbdwp[6] = de1dtdr2;
    dbdwp[7] = de1dtdu2;
    dbdwp[8] = de1dtdp2;

    // Get Transformation Matrix
    dWpdW(dwpdw, W, nx - 2);
    for(int row = 0; row < 3; row++)
    for(int col = 0; col < 3; col++)
        for(int k = 0; k < 3; k++)
        {
            dBodWo[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];
        }
    
    
    // *********************
    // INLET JACOBIANS
    // *********************

    if(u1 < c1)
    {
        i1 = 0;
        i2 = 1;
        r1 = rho[i1];
        r2 = rho[i2];
        p1 = p[i1];
        p2 = p[i2];
        u1 = u[i1];
        u2 = u[i2];
        c1 = c[i1];
        c2 = c[i2];
        t1 = T[i1];

        // Shorthand
        drho = r1 - r2;
        dp = p1 - p2;
        du = u1 - u2;
        cr = r1 * c1;
        uu = u1 * u1;
        gamr = (gam - 1) / (gam + 1);
        fu = 1 - gamr * u1 * u1 / a2;

        // Speed of Sound
        dc1dr1 = - p1 * gam / (2 * cr * r1);
        dc2dr2 = - p2 * gam / (2 * c2 * r2 * r2);
        dc1dp1 = gam / (2 * cr);
        dc2dp2 = gam / (2 * c2 * r2);

        // Eigenvalue
        eig1 = (u1 + u2) / 2;
        eig3 = eig1 - (c1 + c2);

        deig1du1 = 1.0 / 2.0;
        deig1du2 = 1.0 / 2.0;

        deig3dr1 = - dc1dr1 / 2;
        deig3du1 = deig1du1;
        deig3dp1 = - dc1dp1 / 2;
        deig3dr2 = - dc2dr2 / 2;
        deig3du2 = deig1du2;
        deig3dp2 = - dc2dp2 / 2;

        // Riemann Invariants
        R3 = - eig3 * (dp - cr * du);

        dR3dr1 = c1 * du * eig3 + du * eig3 * r1 * dc1dr1 + (-dp + cr * du) * deig3dr1; 
        dR3du1 = - cr * eig3 + (-dp + cr * du) * deig3du1;
        dR3dp1 = - eig3 + du * eig3 * r1 * dc1dp1 + (-dp + cr * du) * deig3dp1;
        dR3dr2 = (-dp + cr * du) * deig3dr2;
        dR3du2 = - cr * eig3 + (-dp + cr * du) * deig3du2;
        dR3dp2 = eig3 + (-dp + cr * du) * deig3dp2;


        // dp1
        double dp1du1, dp1du1du1;
        dp1du1 = -2 * ptin * pow(fu, 1 / (gam - 1)) * gamr * u1 / a2;
        dp1du1du1 = ( 2 * gamr * ptin * pow( fu, ( gam / (gam - 1) - 2) )
                    * ( 1 - gam + gamr * uu * (gam + 1) ) )
                    / ( a2 * (gam - 1) );

        // du1
        du1dt = R3 / (dp1du1 - cr);
        du1dtdr1 = dR3dr1 / (dp1du1 - cr)
                   + R3 * (c1 + r1 * dc1dr1) / pow((dp1du1 - cr), 2);
        du1dtdu1 = dR3du1 / (dp1du1 - cr)
                   - R3 * dp1du1du1 / pow((dp1du1 - cr), 2);
        du1dtdp1 = (R3 * r1 * dc1dp1) / pow((dp1du1 - cr), 2)
                   + dR3dp1 / (dp1du1 - cr);
        du1dtdr2 = dR3dr2 / (dp1du1 - cr);
        du1dtdu2 = dR3du2 / (dp1du1 - cr);
        du1dtdp2 = dR3dp2 / (dp1du1 - cr);

        // dp1/dt
        dp1dt  = dp1du1 * du1dt;
        dp1dtdr1 = dp1du1 * du1dtdr1;
        dp1dtdu1 = du1dt * dp1du1du1 + dp1du1 * du1dtdu1;
        dp1dtdp1 = dp1du1 * du1dtdp1;
        dp1dtdr2 = dp1du1 * du1dtdr2;
        dp1dtdu2 = dp1du1 * du1dtdu2;
        dp1dtdp2 = dp1du1 * du1dtdp2;

        // dt1
        double dt1dp1, dt1dp1dp1;
        dt1dp1 = Ttin / ptin * (gam - 1) / gam * pow(p1 / ptin, - 1 / gam);
        dt1dp1dp1 = - Ttin * (gam - 1) / pow((gam * ptin), 2)
                    * pow(p1 / ptin, - (1 + gam) / gam);

        // dr1/dt
        double dr1dp1, dr1dp1dp1;
        dr1dp1 = 1 / (R * t1) + p1 * dt1dp1 / (R * t1 * t1);
        dr1dp1dp1 = 2 * p1 * dt1dp1 * dt1dp1 / (R * t1 * t1 * t1)
                    - 2 * dt1dp1 / (R * t1 * t1)
                    - p1 * dt1dp1dp1 / (R * t1 * t1);
        dr1dt = dr1dp1 * dp1dt;

        dr1dtdr1 = dr1dp1 * dp1dtdr1;
        dr1dtdu1 = dr1dp1 * dp1dtdu1;
        dr1dtdp1 = dr1dp1 * dp1dtdp1 + dp1dt * dr1dp1dp1;
        dr1dtdr2 = dr1dp1 * dp1dtdr2;
        dr1dtdu2 = dr1dp1 * dp1dtdu2;
        dr1dtdu2 = dr1dp1 * dp1dtdu2;

        // dru1/dt
        // dru1dt = dr1dt * u1 + du1dt * r1;

        dru1dtdr1 = du1dt + u1 * dr1dtdr1 + r1 * du1dtdr1;
        dru1dtdu1 = dr1dt + u1 * dr1dtdu1 + r1 * du1dtdu1;
        dru1dtdp1 = u1 * dr1dtdp1 + r1 * du1dtdp1;
        dru1dtdr2 = u1 * dr1dtdr2 + r1 * du1dtdr2;
        dru1dtdu2 = u1 * dr1dtdu2 + r1 * du1dtdu2;
        dru1dtdu2 = u1 * dr1dtdp2 + r1 * du1dtdp2;

        // de1/dt
        de1dtdr1 = du1dt * u1 + dp1dtdr1 * Cv / R + uu * dr1dtdr1 / 2 + cr * du1dtdr1;
        de1dtdu1 = dr1dt * u1 + du1dt * r1 + dp1dtdu1 * Cv / R
                  + uu * dr1dtdu1 / 2 + cr * du1dtdu1;
        de1dtdp1 = dp1dtdp1 * Cv / R + uu * dr1dtdp1 / 2 + cr * du1dtdp1;
        de1dtdr2 = dp1dtdr2 * Cv / R + u1 + cr * du1dtdr2 + uu * dr1dtdr2 / 2;
        de1dtdu2 = dp1dtdu2 * Cv / R + u1 + cr * du1dtdu2 + uu * dr1dtdu2 / 2;
        de1dtdp2 = dp1dtdp2 * Cv / R + u1 + cr * du1dtdp2 + uu * dr1dtdp2 / 2;


        dbdwp[0] = dr1dtdr1;
        dbdwp[1] = dr1dtdu1;
        dbdwp[2] = dr1dtdp1;
        dbdwp[3] = dru1dtdr1;
        dbdwp[4] = dru1dtdu1;
        dbdwp[5] = dru1dtdp1;
        dbdwp[6] = de1dtdr1;
        dbdwp[7] = de1dtdu1;
        dbdwp[8] = de1dtdp1;

        // Get Transformation Matrix
        dWpdW(dwpdw, W, 0);
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                dBidWi[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];
            }
        
        dbdwp[0] = dr1dtdr2;
        dbdwp[1] = dr1dtdu2;
        dbdwp[2] = dr1dtdp2;
        dbdwp[3] = dru1dtdr2;
        dbdwp[4] = dru1dtdu2;
        dbdwp[5] = dru1dtdp2;
        dbdwp[6] = de1dtdr2;
        dbdwp[7] = de1dtdu2;
        dbdwp[8] = de1dtdp2;

        // Get Transformation Matrix
        dWpdW(dwpdw, W, 1);
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                dBidWd[row * 3 + col] += dbdwp[row * 3 + k] * dwpdw[k * 3 + col];
            }

        }
}
