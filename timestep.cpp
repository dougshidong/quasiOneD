// Time Stepping Schemes

#include<vector>
#include<math.h>
#include<iostream>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include "flux.h"
#include "convert.h"
#include "globals.h"
#include "petscGMRES.h"
#include "residuald1.h"

int ki, kip;

void EulerExplicitStep(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &V,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W);

void rk4(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &V,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W);

void jamesonrk(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &V,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W);

void eulerImplicit(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &V,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W);

void crankNicolson(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &V,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W);

double lineSearchW(
    const std::vector <double> &Resi,
    const std::vector <double> &S,
    const std::vector <double> &W,
    VectorXd &dW);

std::vector <double> V(nx);
std::vector <double> Flux;
std::vector <double> Q;

std::vector <double> W1, W2, W3, Wtemp;
std::vector <double> Resi0, Resi1, Resi2;

void stepInTime(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W)
{
    for(int i = 0; i < nx; i++)
    {
        V[i] = (S[i] + S[i + 1]) / 2.0 * dx[i];
    }
    if(StepScheme == 0)
    {
        EulerExplicitStep(S, dx, V, dt, Resi, W);
    }
    else if(StepScheme == 1)
    {
        rk4(S, dx, V, dt, Resi, W);
    }
    else if(StepScheme == 2)
    {
        jamesonrk(S, dx, V, dt, Resi, W);
    }
    else if(StepScheme == 3)
    {
        eulerImplicit(S, dx, V, dt, Resi, W);
    }
    else if(StepScheme == 4)
    {
        crankNicolson(S, dx, V, dt, Resi, W);
    }
}

void initializeTimeStep(int nx)
{
    V.resize(nx);
    Flux.resize(3 * (nx + 1));
    Q.resize(3 * nx);
    Wtemp.resize(3 * nx);
    Resi1.resize(3 * nx);
    if(StepScheme == 1) // RK4
    {
        W1.resize(3 * nx);
        W2.resize(3 * nx);
        W3.resize(3 * nx);
        Wtemp.resize(3 * nx);
        Resi0.resize(3 * nx);
        Resi1.resize(3 * nx);
        Resi2.resize(3 * nx);
    }
}

// Domain Residual R = FS_i+1/2 - FS_i-1/2 - Qi
void getDomainResi(
    const std::vector <double> &W,
    const std::vector <double> &Flux,
    const std::vector <double> &S,
    std::vector <double> &Resi)
{
    WtoQ(W, Q, S);
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = i * 3 + k;
            kip = (i + 1) * 3 + k;
            Resi[ki] = Flux[kip] * S[i + 1] - Flux[ki] * S[i] - Q[ki];
        }
    }
}

// Euler Explicit
void EulerExplicitStep(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &V,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W)
{
    getFlux(Flux, W);
    getDomainResi(W, Flux, S, Resi);

    for(int k = 0; k < 3; k++)
    for(int i = 1; i < nx - 1; i++)
    {
        ki = i * 3 + k;
        W[ki] = W[ki] - (dt[i] / V[i]) * Resi[ki];
    }

    return;
}

// 4th order Runge - Kutta Stepping Scheme
void rk4(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &V,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W)
{
    // Residual 0
    getFlux(Flux, W);
    getDomainResi(W, Flux, S, Resi0);
    // RK1
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = i * 3 + k;
            W1[ki] = W[ki] - (dt[i] / 2) * Resi0[ki] / dx[i];
        }
        W1[0 * 3 + k] = W[0 * 3 + k];
        W1[(nx - 1) * 3 + k] = W[(nx - 1) * 3 + k];
    }

    // Residual 1
    getFlux(Flux, W1);
    getDomainResi(W1, Flux, S, Resi1);

    // RK2
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = i * 3 + k;
            W2[ki] = W[ki] - (dt[i] / 2) * Resi1[ki] / dx[i];
        }
        W2[0 * 3 + k] = W[0 * 3 + k];
        W2[(nx - 1) * 3 + k] = W[(nx - 1) * 3 + k];
    }

    // Residual 2
    getFlux(Flux, W2);
    getDomainResi(W2, Flux, S, Resi2);

    // RK3
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = i * 3 + k;
            W3[ki] = W[ki] - (dt[i] / 2) * Resi2[ki] / dx[i];
        }
    }

    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = i * 3 + k;
            Wtemp[ki] = ((double)1.0 / 6.0) * (W[ki] + 2 * W1[ki] + 2 * W2[ki] + W3[ki]);
            //Resi[ki] = (2 * Resi0[ki] + 2 * Resi1[ki] + Resi2[ki]) / 6.0;
            Resi[ki] = (Wtemp[ki] - W[ki]) * V[i] / dt[i];
            W[ki] = Wtemp[ki];
        }
    }
}


// Jameson's 4th order Runge - Kutta Stepping Scheme
void jamesonrk(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &V,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W)
{
    // Initialize First Stage
    for(int k = 0; k < 3; k++)
    {
        for(int i = 0; i < nx; i++)
        {
            ki = i * 3 + k;
            Wtemp[ki] = W[ki];
        }
    }
    // 1-4 Stage
    for(int r = 1; r < 5; r++)
    {
        // Get Flux
        getFlux(Flux, Wtemp);
        // Calculate Residuals
        getDomainResi(Wtemp, Flux, S, Resi1);
        // Step in RK time
        for(int k = 0; k < 3; k++)
        {
            for(int i = 1; i < nx - 1; i++)
            {
                ki = i * 3 + k;
                Wtemp[ki] = W[ki] - (dt[i] / (5 - r)) * Resi1[ki] / dx[i];
            }
            Wtemp[0 * 3 + k] = W[0 * 3 + k];
            Wtemp[(nx - 1) * 3 + k] = W[(nx - 1) * 3 + k];
        }
    }

    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = i * 3 + k;
            Resi[ki] = (Wtemp[ki] - W[ki]) * V[i] / dt[i];
            W[ki] = Wtemp[ki];
        }
    }
}

void dFdW(
    std::vector <double> &J,
    double Sp, double Sm,
    double rho, double u, double c)
{
    J[0] = 0.0;
    J[1] = 1.0;
    J[2] = 0.0;
    J[3] = u * u * (gam - 3.0) / 2.0;
    J[4] = u * (3.0 - gam);
    J[5] = gam - 1.0;
    J[6] = ( pow(u, 3) * (gam - 1.0) 
           * (gam - 2.0) - 2.0 * u * c * c ) 
           / (2.0 * (gam - 1.0));
    J[7] = ( 2.0 * c * c + u * u 
           * ( -2.0 * gam * gam + 5.0 * gam - 3.0 ) )
           / (2.0 * (gam - 1.0));
    J[8] = u * gam;

    double Sa = (Sp + Sm) / 2.0;
    for(int i = 0; i < 9; i++)
    {
        J[i] *= Sa;
    }
    
    double dpdw[3];
    dpdw[0] = (gam - 1) / 2.0 * u * u;
    dpdw[1] = - (gam - 1) * u;
    dpdw[2] = (gam - 1);

    J[3] -= dpdw[0] * (Sp - Sm);
    J[4] -= dpdw[1] * (Sp - Sm);
    J[5] -= dpdw[2] * (Sp - Sm);
}

void eulerImplicit(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &V,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W)
{
    // Get Flux
    getFlux(Flux, W);
    // Calculate Residuals
    getDomainResi(W, Flux, S, Resi1);
    Eigen::VectorXd RHS(3 * nx);
    RHS.setZero();
    Eigen::SparseMatrix <double> A(3 * nx, 3 * nx);
    A = evaldRdW(W, dx, dt, S);
//  A = evaldRdW_FD(W, S);
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            A.coeffRef(i,j) = 0.0;
        }
    }
    for(int i = 0; i < 3 * nx; i++)
    {
        A.coeffRef(i,i) += 1.0;
    }

    int Wik;
    for(int Wi = 1; Wi < nx - 1; Wi++)
    {
        for(int Wk = 0; Wk < 3; Wk++)
        {
            Wik = Wi * 3 + Wk;
            RHS[Wik] += - dt[Wi] / V[Wi] * Resi1[Wik];

            if(Wi == 0 || Wi == nx-1)
            {
                RHS[Wik] = -dt[Wi] / V[Wi] * Resi[Wik];
            }
        }
    }
//  std::cout<<A<<std::endl;
//  std::cout<<RHS<<std::endl;
//  *((unsigned int*)0) = 0xDEAD;
    VectorXd Wt(3 * nx);
//  Wt = solveGMRES(A, RHS);
    SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver1;
    slusolver1.compute(A);
    if(slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
    Wt = slusolver1.solve(RHS);
    double currentR = 0;
    for(int i = 0; i < nx; i++)
        currentR += Resi[i * 3 + 0] * Resi[i * 3 + 0];
    currentR = sqrt(currentR);
    double alpha = 1;
    if(1.0/currentR > 1e2)
    {
        alpha = 1.25 * pow(log10(1.0/currentR/1e1), 3);
        std::cout<<alpha<<std::endl;
    }
    std::cout<<alpha<<std::endl;
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = i * 3 + k;
            W[ki] += Wt[ki] * alpha;
            Resi[ki] = Wt[ki] * V[i] / dt[i];
        }
    }
}
void crankNicolson(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &V,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W)
{
    // Get Flux
    getFlux(Flux, W);
    // Calculate Residuals
    getDomainResi(W, Flux, S, Resi1);
    Eigen::VectorXd RHS(3 * nx);
    RHS.setZero();
    Eigen::SparseMatrix <double> A(3 * nx, 3 * nx);
    A = 0.5 * evaldRdW(W, dx, dt, S);
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            A.coeffRef(i,j) = 0.0;
        }
    }
    for(int i = 0; i < 3 * nx; i++)
    {
        A.coeffRef(i,i) += 1.0;
    }

    int Wik;
    for(int Wi = 1; Wi < nx - 1; Wi++)
    {
        for(int Wk = 0; Wk < 3; Wk++)
        {
            Wik = Wi * 3 + Wk;
            RHS[Wik] += - dt[Wi] / V[Wi] * Resi1[Wik];

            if(Wi == 0 || Wi == nx-1)
            {
                RHS[Wik] = -dt[Wi] / V[Wi] * Resi[Wik];
            }
        }
    }
//  std::cout<<A<<std::endl;
//  std::cout<<RHS<<std::endl;
//  *((unsigned int*)0) = 0xDEAD;
    VectorXd Wt(3 * nx);
//  Wt = solveGMRES(A, RHS);
    SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver1;
    slusolver1.compute(A);
    if(slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
    Wt = slusolver1.solve(RHS);
    double currentR = 0;
    for(int i = 0; i < nx; i++)
        currentR += Resi[i * 3 + 0] * Resi[i * 3 + 0];
    currentR = sqrt(currentR);
    double alpha = 1;
    if(1.0/currentR > 1e2)
    {
        alpha = 1.1 * pow(log10(1.0/currentR/1e1), 3.0);
        std::cout<<alpha<<std::endl;
    }
    std::cout<<alpha<<std::endl;
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = i * 3 + k;
            W[ki] += Wt[ki] * alpha;
            Resi[ki] = Wt[ki] * V[i] / dt[i];
        }
    }
}

double lineSearchW(
    const std::vector <double> &Resi,
    const std::vector <double> &S,
    const std::vector <double> &W,
    VectorXd &dW)
{
    double alpha = 1;
    double c1 = 1e-4;
    std::vector <double> tempS(nx + 1);

    double c_pk_grad = 0;
    c_pk_grad = c1 * dW.dot(dW);

    double currentVal = 0;
    for(int i = 0; i < 3 * nx; i++)
        currentVal += Resi[i] * Resi[i];
    currentVal = sqrt(currentVal);

    std::vector <double> tempD(nDesVar);
    for(int i = 0; i < 3 * nx; i++)
        Wtemp[i] = W[i] + alpha * dW[i];
    double newVal = 0;
    getFlux(Flux, Wtemp);
    getDomainResi(Wtemp, Flux, S, Resi1);
    for(int i = 0; i < 3 * nx; i++)
        newVal += Resi1[i] * Resi1[i];
    newVal = sqrt(newVal);


    std::cout<<newVal<<std::endl;
    std::cout<<currentVal<<std::endl;
    std::cout<<currentVal + alpha * c_pk_grad<<std::endl;
    while(newVal > (currentVal + alpha * c_pk_grad) && alpha > 1e-16)
    {
        alpha = alpha * 0.75;
//      std::cout<<"Alpha Reduction: "<<alpha<<std::endl;

        for(int i = 0; i < 3 * nx; i++)
            Wtemp[i] = W[i] + alpha * dW[i];
        getFlux(Flux, Wtemp);
        getDomainResi(Wtemp, Flux, S, Resi1);
        newVal = 0;
        for(int i = 0; i < 3 * nx; i++)
            newVal += Resi1[i] * Resi1[i];
        newVal = sqrt(newVal);

//      std::cout<<"newVal: "<<newVal<<std::endl;
//      std::cout<<"currentVal + alpha/2.0 * c_pk_grad: "<<
//      currentVal + alpha/ 2.0 * c_pk_grad<<std::endl;
    }
//  if(alpha < 1e-16) std::cout<<"Error. Can't find step size"<<std::endl;

    std::cout<<"Alpha: "<<alpha<<std::endl;

    return alpha;
}
