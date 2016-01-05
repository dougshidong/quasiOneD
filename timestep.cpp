// Time Stepping Schemes

#include<vector>
#include<math.h>
#include "flux.h"
#include "globals.h"
#include "flovar.h"

// Euler Explicit
void EulerExplicitStep(std::vector <double> S,
                       std::vector <double> V,
                       std::vector <double> dt,
                       std::vector <double> Flux,
                       std::vector <double> Q,
                       std::vector <double> &Resi,
                       std::vector <double> &W)
{
    int ki;
    for(int k = 0; k < 3; k++)
    {
        for(int i = 1; i < nx - 1; i++)
        {
            ki = k * nx + i;
//            std::cout<<Flux[ki + 1]<<std::endl;
  //          std::cout<<Flux[ki]<<std::endl;
            Resi[ki] = Flux[ki + 1] * S[i + 1]
                 - Flux[ki] * S[i] - Q[ki];
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


// Jameson's 4th order Runge - Kutta Stepping Scheme
void rk4(std::vector <double> dx, std::vector <double> S, 
         std::vector <double> dt, 
         std::vector <double> &W,
         std::vector <double> Q,    
         std::vector <double> &Resi,
         std::vector <double> &Flux)
{
    double ki;
    std::vector <double> Resi0(3 * nx, 0), Resi1(3 * nx, 0), Resi2(3 * nx, 0);
    std::vector <double> W1(3 * nx, 0), W2(3 * nx, 0), W3(3 * nx, 0);
    std::vector <double> Q1(3 * nx, 0), Q2(3 * nx, 0);
    std::vector <double> utemp(nx), rhotemp(nx), ptemp(nx), ctemp(nx);

    std::vector <double> ResiR(3 * nx * 3, 0);
    std::vector <double> WR(3 * nx * 4, 0);

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

    Flux_StegerWarmingV(Flux, W1, utemp, ctemp, rhotemp);

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

    Flux_StegerWarmingV(Flux, W2, utemp, ctemp, rhotemp);

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

    return;
}
