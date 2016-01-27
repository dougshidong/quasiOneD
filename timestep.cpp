// Time Stepping Schemes

#include<vector>
#include<math.h>
#include "flux.h"
#include "convert.h"
#include "globals.h"
#include "flovar.h"

std::vector <double> Flux(3 * (nx + 1), 0);
// Euler Explicit
void EulerExplicitStep(std::vector <double> S,
                       std::vector <double> V,
                       std::vector <double> dt,
                       std::vector <double> Q,
                       std::vector <double> &Resi,
                       std::vector <double> &W,
                       std::vector <double> F)
{
    int ki;
    getFlux(Flux, W, F);
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
            Resi[ki] = Flux[ki + 1] * S[i + 1] - Flux[ki] * S[i] - Q[ki];
        }
        Resi[k * nx + 0] = 0;
        Resi[k * nx + (nx - 1)] = 0;
    }
    for(int k = 0; k < 3; k++)
    for(int i = 1; i < nx - 1; i++)
    {
        ki = k * nx + i;
        W[ki] = W[ki] - (dt[i] / V[i]) * Resi[ki];
    }

    return;
}

std::vector <double> Resi0(3 * nx, 0), Resi1(3 * nx, 0), Resi2(3 * nx, 0);
std::vector <double> W1(3 * nx, 0), W2(3 * nx, 0), W3(3 * nx, 0);
std::vector <double> F1(3 * nx, 0), F2(3 * nx, 0);
std::vector <double> Q1(3 * nx, 0), Q2(3 * nx, 0);
std::vector <double> utemp(nx), rhotemp(nx), ptemp(nx), ctemp(nx);
// 4th order Runge - Kutta Stepping Scheme
void rk4(std::vector <double> dx, 
         std::vector <double> S, 
         std::vector <double> dt, 
         std::vector <double> &W,
         std::vector <double> F,
         std::vector <double> Q,    
         std::vector <double> &Resi)
{
    double ki;

    getFlux(Flux, W, F);
    // Residual 0
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
            Resi0[ki] = Flux[ki + 1] * S[i + 1] - Flux[ki] * S[i] - Q[ki];
        }
    }
    // RK1
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
            W1[ki] = W[ki] - (dt[i] / 2) * Resi0[ki] / dx[i];
        }
        W1[k * nx] = W[k * nx + 0];
        W1[k * nx + (nx - 1)] = W[k * nx + (nx - 1)];
    }
    for(int i = 0; i < nx; i++)
    {
        rhotemp[i] = W1[i];
        utemp[i] = W1[1 * nx + i] / rhotemp[i];
        ptemp[i] = (gam - 1) * (W1[2 * nx + i] - rhotemp[i] * utemp[i] / 2);
        ctemp[i] = sqrt(gam * ptemp[i] / rhotemp[i]);

        Q1[nx + i] = ptemp[i] * (S[i + 1] - S[i]);
    }

    WtoF(W1, F1);
    getFlux(Flux, W1, F1);

    // Residual 1
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
            Resi1[ki] = Flux[ki + 1] * S[i + 1] - Flux[ki] * S[i] - Q1[ki];
        }
    }

    // RK2
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
            W2[ki] = W[ki] - (dt[i] / 2) * Resi1[ki] / dx[i];
        }
        W2[k * nx] = W[k * nx + 0];
        W2[k * nx + (nx - 1)] = W[k * nx + (nx - 1)];
    }
    for(int i = 0; i < nx; i++)
    {
        rhotemp[i] = W1[i];
        utemp[i] = W1[1 * nx + i] / rhotemp[i];
        ptemp[i] = (gam - 1) * (W1[2 * nx + i] - rhotemp[i] * utemp[i] / 2);
        ctemp[i] = sqrt(gam * ptemp[i] / rhotemp[i]);

        Q2[nx + i] = ptemp[i] * (S[i + 1] - S[i]);
    }

    WtoF(W2, F2);
    getFlux(Flux, W2, F2);

    // Residual 2
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
            Resi2[ki] = Flux[ki + 1] * S[i + 1] - Flux[ki] * S[i] - Q2[ki];
        }
    }

    // RK3
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
            W3[ki] = W[ki] - (dt[i] / 2) * Resi2[ki] / dx[i];
        }
    }

    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
            W[ki] = ((double)1.0 / 6.0) * (W[ki] + 2 * W1[ki] + 2 * W2[ki] + W3[ki]);
            Resi[ki] = (2 * Resi0[ki] + 2 * Resi1[ki] + Resi2[ki]) / 6.0;
        }
    }
}


std::vector <double> Resitemp(3 * nx, 0);
std::vector <double> Wtemp(3 * nx, 0);
std::vector <double> Ftemp(3 * nx, 0);
std::vector <double> Qtemp(3 * nx, 0);
// Jameson's 4th order Runge - Kutta Stepping Scheme
void jamesonrk(std::vector <double> dx, 
         std::vector <double> S,
         std::vector <double> V,
         std::vector <double> dt, 
         std::vector <double> &W,
         std::vector <double> F,
         std::vector <double> &Resi)
{
    double ki;

    // Initialize First Stage
    for(int k = 0; k < 3; k++)
    {
        for(int i = 0; i < nx; i++)
        {
            ki = k * nx + i;
            Wtemp[ki] = W[ki];
            Ftemp[ki] = F[ki];
        }
    }
    for(int i = 0; i < nx; i++)
    {
        rhotemp[i] = W[i];
        utemp[i] = W[1 * nx + i] / rhotemp[i];
        ptemp[i] = (gam - 1) * (W[2 * nx + i] - rhotemp[i] * utemp[i] / 2);
        ctemp[i] = sqrt(gam * ptemp[i] / rhotemp[i]);

        Qtemp[nx + i] = ptemp[i] * (S[i + 1] - S[i]);
    }
    // 1-4 Stage
    for(int r = 1; r < 5; r++)
    {
        // Get Flux
        getFlux(Flux, Wtemp, Ftemp);
        // Calculate Residuals
        for(int k = 0; k < 3; k++)
        {
            for(int i = 1; i < nx - 1; i++)
            {
                ki = k * nx + i;
                Resitemp[ki] = Flux[ki + 1] * S[i + 1] - Flux[ki] * S[i] - Qtemp[ki];
            }
        }
        // Step in RK time
        for(int k = 0; k < 3; k++)
        {
            for(int i = 1; i < nx - 1; i++)
            {
                ki = k * nx + i;
                Wtemp[ki] = W[ki] - (dt[i] / (5 - r)) * Resitemp[ki] / dx[i];
            }
            Wtemp[k * nx + 0] = W[k * nx + 0];
            Wtemp[k * nx + (nx - 1)] = W[k * nx + (nx - 1)];
        }
        // Calculate temporary variables
        for(int i = 0; i < nx; i++)
        {
            rhotemp[i] = Wtemp[i];
            utemp[i] = Wtemp[1 * nx + i] / rhotemp[i];
            ptemp[i] = (gam - 1) * (Wtemp[2 * nx + i] - rhotemp[i] * utemp[i] / 2);
            ctemp[i] = sqrt(gam * ptemp[i] / rhotemp[i]);

            Qtemp[nx + i] = ptemp[i] * (S[i + 1] - S[i]);
        }
    }

    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
            Resi[ki] = (Wtemp[ki] - W[ki]) * dt[i] / V[i];
            W[ki] = Wtemp[ki];
        }
    }
}
