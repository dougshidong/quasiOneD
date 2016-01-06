// Convert W to F or primitive variables

#include<vector>
#include<math.h>
#include"globals.h"
#include"flovar.h"


// Get primitive variables
void WtoP(std::vector <double> W,
          std::vector <double> &rho,
          std::vector <double> &u,
          std::vector <double> &e)
{
    for(int i = 0; i < nx; i++)
    {
        rho[i] = W[0 * nx + i];
        u[i] = W[1 * nx + i] / rho[i];
        e[i] = W[2 * nx + i];
    }
}

// Get more primitive variables
void WtoP(std::vector <double> W,
          std::vector <double> &rho,
          std::vector <double> &u,
          std::vector <double> &e,
          std::vector <double> &p,
          std::vector <double> &c)
{
    for(int i = 0; i < nx; i++)
    {
        rho[i] = W[0 * nx + i];
        u[i] = W[1 * nx + i] / rho[i];
        e[i] = W[2 * nx + i];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2 );
        c[i] = sqrt( gam * p[i] / rho[i] );
    }
}
// Get F
void WtoF(std::vector <double> W,
          std::vector <double> &F)
{
    double w1, w2, w3;
    for(int i = 0; i < nx; i++)
    {
        w1 = W[0 * nx + i];
        w2 = W[1 * nx + i];
        w3 = W[2 * nx + i];
        F[0 * nx + i] = w2;
        F[1 * nx + i] = w2 * w2 / w1 + (gam - 1) * ( w3 - w2 * w2 / (2 * w1) );
        F[2 * nx + i] = ( w3 + (gam - 1) * (w3 - w2 * w2 / (2 * w1)) ) * w2 / w1;
    }
}
