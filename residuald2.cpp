#include<Eigen/Core>
#include<Eigen/Sparse>
#include<math.h>
#include<iostream>
#include"globals.h"
#include"convert.h"
#include"flux.h"
#include"quasiOneD.h"
#include"residuald1.h"

using namespace Eigen;
  
std::vector <MatrixXd> evalddFdWdW(double rho, double u, double p);

void evalddFluxdWdW(
    std::vector <MatrixXd> &ddFluxdWdW1,
    std::vector <MatrixXd> &ddFluxdWdW2,
    std::vector <MatrixXd> &ddFluxdWdW3,
    std::vector <double> W);

// Calculates Residual Hessian
std::vector < SparseMatrix <double> > evalddRdWdW(
    std::vector <double> W,
    std::vector <double> S)
{
    std::vector <SparseMatrix <double> > ddRdWdW(3 * nx);
    SparseMatrix <double> ddRdWd_Rik(3 * nx, 3 * nx);


    // Flux Hessian 
    // ddFluxdWdW1 = Flux(i+1/2) / W(i, i)
    // ddFluxdWdW2 = Flux(i+1/2) / W(i + 1, i + 1)
    // ddFluxdWdW3 = Flux(i+1/2) / W(i, i + 1)
    // ddFluxdWdW4 = ddFluxdWdW3 (not required)
    std::vector <MatrixXd> 
        ddFluxdWdW1(3 * nx), ddFluxdWdW2(3 * nx), ddFluxdWdW3(3 * nx);

    evalddFluxdWdW(ddFluxdWdW1, ddFluxdWdW2, ddFluxdWdW3, W);
    int Rki, Rkip;
    int Wi0, Wi1, Wi2;
    SparseMatrix <double> ddRdWdW_Rik(3 * nx, 3 * nx);
    for(int Ri = 1; Ri < nx - 1; Ri++)
    {
        for(int Rk = 0; Rk < 3; Rk++)
        {
            Rki = Ri * 3 + Rk;
            Rkip = (Ri + 1) * 3 + Rk;

            // R(i) = Flux(i) * S(i) + Flux(i-1) * S(i-1) - Q(i)

            //ddR(i) / dW(i-1)dW(i-1)
            // = ddFlux(i-1) / dW(i-1)dW(i-1) * S(i-1)
            ddRdWdW_Rik.setZero();
            for(int row; row < 3; row++)
            for(int col; col < 3; col++)
            {
                ddRdWdW_Rik.insert(row, col) = -ddFluxdWdW1[Rki].coeffRef(row, col) * S[Ri];
            }
            ddRdWdW[Rki] = ddRdWdW_Rik;

            //ddR(i) / dW(i+1)dW(i+1)
            // = ddFlux(i+1) / dW(i+1)dW(i+1) * S(i+1)
            ddRdWdW_Rik.setZero();
            for(int row; row < 3; row++)
            for(int col; col < 3; col++)
            {
                ddRdWdW_Rik.insert(row, col) = -ddFluxdWdW1[Rkip].coeffRef(row, col) * S[Ri + 1];
            }
            ddRdWdW[Rki] = ddRdWdW_Rik;

            //ddR(i) / dW(i)dW(i+1)
            Wi1 = Ri;
            Wi2 = Ri + 1;

            //ddR(i) / dW(i+1)dW(i)
            Wi1 = Ri + 1;
            Wi2 = Ri;
        }

    }
    return ddRdWdW;
}

std::vector <SparseMatrix <double> > evalddRdWdS_FD(
    std::vector <double> W,
    std::vector <double> S)
{
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);
    SparseMatrix <double> dRdWdSi(3 * nx, nx + 1);

    SparseMatrix <double> ddRidWdS(3 * nx, nx + 1);
    for(int Ri = 0; Ri < 3 * nx; Ri++)
        ddRdWdS[Ri] = ddRidWdS;

    double Resi0, Resi1, Resi2, Resi3, Resi4;
    std::vector <double> Flux(3 * (nx + 1), 0);
    std::vector <double> Wd(3 * nx, 0);
    std::vector <double> Sd(nx + 1, 0);
    std::vector <double> Q(3 * nx, 0);
    double h = 0.0000001;
    double pertW, pertS;
    int ki, kip;
    for(int Ri = 1; Ri < nx - 1; Ri++)
    {
        std::cout<<"Ri is "<<Ri<<std::endl;
        for(int k = 0; k < 3; k++)
        {
            ki = Ri * 3 + k;
            kip = (Ri + 1) * 3 + k;
            
            WtoQ(W, Q, S);
            getFlux(Flux, W);
            Resi0 = Flux[kip] * S[Ri + 1] - Flux[ki] * S[Ri] - Q[ki];

            for(int Wi = 0; Wi < 3 * nx; Wi++)
            {
                pertW = W[Wi] * h;
                for(int Si = 0; Si < nx + 1; Si++)
                {
                    for(int m = 0; m < 3 * nx; m++)
                        Wd[m] = W[m];
                    for(int m = 0; m < nx + 1; m++)
                        Sd[m] = S[m];
                    // R1
                    Wd[Wi] = W[Wi];
                    Sd[Si] = S[Si];

                    pertS = S[Si] * h;
                    Wd[Wi] = W[Wi] + pertW;
                    Sd[Si] = S[Si] + pertS;

                    WtoQ(Wd, Q, Sd);
                    getFlux(Flux, Wd);
                    
                    Resi1 = Flux[kip] * Sd[Ri + 1] - Flux[ki] * Sd[Ri] - Q[ki];
                    
                    // R2
                    Wd[Wi] = W[Wi];
                    Sd[Si] = S[Si];

                    pertS = S[Si] * h;
                    Wd[Wi] = W[Wi] + pertW;
                    Sd[Si] = S[Si] - pertS;

                    WtoQ(Wd, Q, Sd);
                    getFlux(Flux, Wd);

                    Resi2 = Flux[kip] * Sd[Ri + 1] - Flux[ki] * Sd[Ri] - Q[ki];
                    
                    // R3
                    Wd[Wi] = W[Wi];
                    Sd[Si] = S[Si];

                    pertS = S[Si] * h;
                    Wd[Wi] = W[Wi] + pertW;
                    Sd[Si] = S[Si] - pertS;

                    WtoQ(Wd, Q, Sd);
                    getFlux(Flux, Wd);

                    Resi3 = Flux[kip] * Sd[Ri + 1] - Flux[ki] * Sd[Ri] - Q[ki];

                    // R4
                    Wd[Wi] = W[Wi];
                    Sd[Si] = S[Si];

                    pertS = S[Si] * h;
                    Wd[Wi] = W[Wi] + pertW;
                    Sd[Si] = S[Si] - pertS;

                    WtoQ(Wd, Q, Sd);
                    getFlux(Flux, Wd);

                    Resi4 = Flux[kip] * Sd[Ri + 1] - Flux[ki] * Sd[Ri] - Q[ki];

                    ddRdWdS[ki].insert(Wi, Si) = (Resi1 - Resi2 - Resi3 + Resi4)
                                                 / (4 * pertW * pertS);
                }// Si Loop
            }// Wi Loop
        }// k Loop
    }// Ri Loop
    return ddRdWdS;
}

std::vector <SparseMatrix <double> > evalddRdWdS(
    std::vector <double> W,
    std::vector <double> S)
{
    // Allocate ddRdWdS Sparse Matrix
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);
    SparseMatrix <double> dRdWdSi(3 * nx, nx + 1);
    dRdWdSi.reserve(3 * 9 * (nx - 2));
    for(int Si = 0; Si < nx + 1; Si++)
        ddRdWdS[Si] = dRdWdSi;

    // Get Jacobians and Fluxes
    std::vector <double> Ap(nx * 3 * 3, 0), An(nx * 3 * 3, 0);
    if(FluxScheme == 1) ScalarJac(W, Ap, An);
    // Evaluate dpdW
    std::vector <double> dpdW(3 * nx, 0);
    dpdW = evaldpdW(W, S);

    int Wi, k, rowi, coli;
    double val;
    for(int Si = 1; Si < nx; Si++)
    {
        for(int Ri = 1; Ri < nx - 1; Ri++)
        {
            Wi = Ri - 1;
            if(Ri == Si && Wi >= 0)
            {
                for(int row = 0; row < 3; row++)
                for(int col = 0; col < 3; col++)
                {
                    k = row * 3 + col;
                    rowi = Ri * 3 + row;
                    coli = Wi * 3 + col;

                    val = - Ap[Wi * 9 + k];
                    ddRdWdS[Si].insert(rowi, coli) = val;
                }
            }

            Wi = Ri;
            if(Ri == Si && Wi >= 0 && Wi <= nx - 1)
            {
                for(int row = 0; row < 3; row++)
                for(int col = 0; col < 3; col++)
                {
                    k = row * 3 + col;
                    rowi = Ri * 3 + row;
                    coli = Wi * 3 + col;

                    val -= An[Wi * 9 + k];
                    if(row == 1) 
                    {
                        val += dpdW[Wi * 3 + col];
                    }

                    ddRdWdS[Si].insert(rowi, coli) = val;
                }
            }

            if((Ri + 1) == Si && Wi >= 0 && Wi <= nx - 1)
            {
                for(int row = 0; row < 3; row++)
                for(int col = 0; col < 3; col++)
                {
                    k = row * 3 + col;
                    rowi = Ri * 3 + row;
                    coli = Wi * 3 + col;

                    val = Ap[Wi * 9 + k];
                    if(row == 1) 
                    {
                        val -= dpdW[Wi * 3 + col];
                    }

                    ddRdWdS[Si].insert(rowi, coli) = val;
                }
            }

            Wi = Ri + 1;
            if(Ri + 1 == Si && Wi <= nx - 1)
            {
                for(int row = 0; row < 3; row++)
                for(int col = 0; col < 3; col++)
                {
                    k = row * 3 + col;
                    rowi = Ri * 3 + row;
                    coli = Wi * 3 + col;

                    val = An[Wi * 9 + k];
                    ddRdWdS[Si].insert(rowi, coli) = val;
                }
            }
        } // Ri Loop
    } // Si Loop
    return ddRdWdS;
}

std::vector <MatrixXd> evalddQdWdW(std::vector <double> W)
{
    std::vector <MatrixXd> ddQdWdW(3 * nx);
    MatrixXd ddQidWdW(3, 3);
    ddQidWdW.setZero();
    for(int i = 0; i < nx; i++)
    {
        ddQdWdW[i * 3 + 0] = ddQidWdW;
        ddQdWdW[i * 3 + 2] = ddQidWdW;
    }
    double rho, u;
    for(int i = 0; i < nx; i++)
    {
        rho = W[i * 3 + 0];
        u = W[i * 3 + 1] / rho;
        ddQidWdW(0, 0) = u * u * (1.0 - gam) / rho;
        ddQidWdW(1, 1) = (1.0 - gam) / rho;
        ddQidWdW(0, 1) = u * (gam - 1.0) / rho;
        ddQidWdW(1, 0) = ddQidWdW(0, 1);
    }
    return ddQdWdW;
}
std::vector <MatrixXd> evalddFdWdW(double rho, double u, double p)
{
    std::vector <MatrixXd> ddFdWdW(3);
    MatrixXd ddfidWdW(3, 3);
    // F1
    ddfidWdW.setZero();
    ddFdWdW[0] = ddfidWdW;
    // F2
    ddfidWdW(0, 0) = -(u * u * (gam - 3.0) / rho);
    ddfidWdW(1, 1) = (3.0 - gam) / rho;
    ddfidWdW(0, 1) = (u * u * (gam - 3.0) / rho);
    ddfidWdW(1, 0) = ddfidWdW(0, 1);
    ddFdWdW[1] = ddfidWdW;
    // F3
    ddfidWdW(0, 0) = ( 2.0 * p * u * gam 
                    + rho * pow(u, 3.0) 
                    * (- 2.0 * gam * gam + 5.0 * gam - 3.0) )
                    / (pow(rho, 2.0) * (gam - 1.0));
    ddfidWdW(1, 1) = 3.0 * u * (1.0 - gam) / rho;
    ddfidWdW(2, 2) = 0;

    ddfidWdW(0, 1) = (- 2.0 * p * gam 
                    + u * u * (5.0 * gam * gam - 11.0 * gam + 6.0) )
                    / (2.0 * pow(rho, 2.0) * (gam - 1.0));
    ddfidWdW(1, 0) = ddfidWdW(0, 1);

    ddfidWdW(0, 2) = -u * gam / rho;
    ddfidWdW(2, 0) = ddfidWdW(0, 2);

    ddfidWdW(1, 2) = gam / rho;
    ddfidWdW(2, 1) = ddfidWdW(0, 2);

    ddFdWdW[2] = ddfidWdW;

    return ddFdWdW;
}

MatrixXd evalddlambdadWdW(std::vector <double> W, int i)
{
    double rho, u, p, c;
    rho = W[i * 3 + 0];
    u = W[i * 3 + 1] / rho;
    p = (gam - 1) * ( W[i * 3 + 2] - rho * u * u / 2.0 );
    c = sqrt(p * gam / rho);

    MatrixXd ddlambdadWdWp(3, 3);

    ddlambdadWdWp(0, 0) = ( 6.0 * p * gam
                          + rho * u * (8.0 * c - u * gam * (gam - 1.0)) )
                          / (16.0 * c * pow(rho, 3.0));
    ddlambdadWdWp(0, 1) = 0.25 * (c * u * (gam - 1.0) / p - 2.0 / rho);
    ddlambdadWdWp(0, 2) = - ( gam * gam * ( 2.0 * p + u * u * (gam - 1.0) * rho ) )
                          / ( 16.0 * pow(c * rho, 3) );

    ddlambdadWdWp(1, 0) = (-4.0 + u * gam * (gam - 1) / c)
                          / (8.0 * rho * rho);
    ddlambdadWdWp(1, 1) = -c * (gam - 1.0) / (4.0 * p);
    ddlambdadWdWp(1, 2) = c * u * (gam - 1.0) / (8.0 * p * p);

    ddlambdadWdWp(2, 0) = - gam * (gam - 1.0)
                          / (8.0 * c * rho * rho);
    ddlambdadWdWp(2, 1) = 0.0;
    ddlambdadWdWp(2, 2) = c * (1.0 - gam) / (8.0 * p * p); 

    return ddlambdadWdWp * dWpdW(W, i);
}

void evalddScalarFluxdWdW(
    std::vector <MatrixXd> &ddFluxdWdW1,
    std::vector <MatrixXd> &ddFluxdWdW2,
    std::vector <MatrixXd> &ddFluxdWdW3,
    std::vector <double> W)
{
    std::vector <double> rho(nx), u(nx), e(nx);
    std::vector <double> T(nx), p(nx), c(nx), Mach(nx);
    WtoP(W, rho, u, e, p, c, T); 
    
    std::vector <MatrixXd> ddFdWdW(3 * nx);
    std::vector <MatrixXd> ddFidWdW(3);

    std::vector <Vector3d> dlambdadW(nx);

    std::vector <MatrixXd> ddlambdadWdW(nx);
    for(int Wi = 0; Wi < nx; Wi++)
    {
        // Evaluate Convective Hessian
        ddFidWdW = evalddFdWdW(rho[Wi], u[Wi], p[Wi]);
        for(int k = 0; k < 3; k++)
        {
            ddFdWdW[Wi * 3 + k] = ddFidWdW[k];
        }
        // Evaluate dlambdadW 
        dlambdadW[Wi] = Map<Vector3d>(evaldlambdadW(W, Wi).data());

        // Evaluate ddlambdadWdW
        ddlambdadWdW[Wi] = evalddlambdadWdW(W, Wi);
    }

    int Fluxkim, Fluxki, Fluxkip;
    MatrixXd ddFluxidWdW(3, 3);
    for(int Fluxi = 1; Fluxi < nx; Fluxi++)
    {
        for(int Fluxk = 0; Fluxk < 3; Fluxk++)
        {
            Fluxkim = (Fluxi - 1) * 3 + Fluxk;
            Fluxki = Fluxi * 3 + Fluxk;
            Fluxkip = (Fluxi + 1) * 3 + Fluxk;
            // ddFluxdWdW1 = F(i+1/2) / W(i, i)
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddFluxidWdW(row, col) = 0.5 * ddFdWdW[Fluxkim](row, col)
                    - 0.5 * Scalareps * ddlambdadWdW[Fluxi - 1](row, col) 
                    * (W[Fluxki] - W[Fluxkim]);
                if(row == col)
                {
                    ddFluxidWdW(row, col) += dlambdadW[Fluxi - 1](row);
                }
            }
            ddFluxdWdW1[Fluxki] = ddFluxidWdW;

            // ddFluxdWdW2 = F(i+1/2) / W(i + 1, i + 1)
            for(int row = 0; row < 3; row++)
            for(int col = 0; col < 3; col++)
            {
                ddFluxidWdW(row, col) = 0.5 * ddFdWdW[Fluxki](row, col)
                    - 0.5 * Scalareps * ddlambdadWdW[Fluxi - 1](row, col) 
                    * (W[Fluxki] - W[Fluxkim]);
                if(row == col)
                {
                    ddFluxidWdW(row, col) -= dlambdadW[Fluxi](row);
                }
            }
            ddFluxdWdW2[Fluxki] = ddFluxidWdW;

            // ddFluxdWdW3 = F(i+1/2) / W(i, i + 1)
            ddFluxidWdW.setZero();
            for(int row = 0; row < 3; row++)
            {
                ddFluxidWdW(row, row) = dlambdadW[Fluxi](row) - dlambdadW[Fluxi - 1](row);
            }
            ddFluxdWdW3[Fluxki] = ddFluxidWdW;
        }
    }
    
}


void evalddFluxdWdW(
    std::vector <MatrixXd> &ddFluxdWdW1,
    std::vector <MatrixXd> &ddFluxdWdW2,
    std::vector <MatrixXd> &ddFluxdWdW3,
    std::vector <double> W)
{
    if(FluxScheme == 1)
    {
        evalddScalarFluxdWdW(ddFluxdWdW1, ddFluxdWdW2, ddFluxdWdW3, W);
    }
}
