// Time Stepping Schemes

#include<vector>
#include<math.h>
#include<iostream>
#include "flux.h"
#include "convert.h"
#include "globals.h"

int ki, kip;
std::vector <double> Q(3 * nx);
std::vector <double> Flux(3 * (nx + 1), 0);
std::vector <double> V(nx);
std::vector <double> Resi0(3 * nx, 0), Resi1(3 * nx, 0), Resi2(3 * nx, 0);
std::vector <double> W1(3 * nx, 0), W2(3 * nx, 0), W3(3 * nx, 0), Wtemp(3 * nx, 0);

void EulerExplicitStep(
    std::vector <double> S,
    std::vector <double> dx,
    std::vector <double> dt,
    std::vector <double> &Resi,
    std::vector <double> &W);

void rk4(
    std::vector <double> S,
    std::vector <double> dx,
    std::vector <double> dt,
    std::vector <double> &Resi,
    std::vector <double> &W);

void jamesonrk(
    std::vector <double> S,
    std::vector <double> dx,
    std::vector <double> dt,
    std::vector <double> &Resi,
    std::vector <double> &W);

void stepInTime(
    std::vector <double> S,
    std::vector <double> dx,
    std::vector <double> dt,
    std::vector <double> &Resi,
    std::vector <double> &W)
{
    for(int i = 0; i < nx; i++)
        V[i] = (S[i] + S[i + 1]) / 2 * dx[i];
    if(StepScheme == 0)
    {
        EulerExplicitStep(S, dx, dt, Resi, W);
    }
    else if(StepScheme == 1)
    {
        rk4(S, dx, dt, Resi, W);
    }
    else if(StepScheme == 2)
    {
        jamesonrk(S, dx, dt, Resi, W);
    }
}

// Domain Residual R = FS_i+1/2 - FS_i-1/2 - Qi
void getDomainResi(
    std::vector <double> W,
    std::vector <double> Flux,
    std::vector <double> S,
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
    std::vector <double> S,
    std::vector <double> dx,
    std::vector <double> dt,
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
    std::vector <double> S,
    std::vector <double> dx,
    std::vector <double> dt,
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
    std::vector <double> S,
    std::vector <double> dx,
    std::vector <double> dt,
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

