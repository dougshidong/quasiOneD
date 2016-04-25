// Time Stepping Schemes

#include<vector>
#include<math.h>
#include<iostream>
#include "flux.h"
#include "convert.h"
#include "globals.h"

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
}

void initializeTimeStep(int nx)
{
    V.resize(nx);
    Flux.resize(3 * (nx + 1));
    Q.resize(3 * nx);
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
    else if(StepScheme == 2) // JRK4
    {
        Wtemp.resize(3 * nx);
        Resi1.resize(3 * nx);
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
    std::vector <double> Flux(3 * (nx + 1));
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
            Resi[ki] = (Wtemp[ki] - W[ki]) * dt[i] / V[i];
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
            Resi[ki] = (Wtemp[ki] - W[ki]) * dt[i] / V[i];
            W[ki] = Wtemp[ki];
        }
    }
}

