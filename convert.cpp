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
          std::vector <double> &c,
          std::vector <double> &T)
{
    for(int i = 0; i < nx; i++)
    {
        rho[i] = W[0 * nx + i];
        u[i] = W[1 * nx + i] / rho[i];
        e[i] = W[2 * nx + i];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2 );
        c[i] = sqrt( gam * p[i] / rho[i] );
        T[i] = p[i] / (rho[i] * R);
    }
}
// Inverse or dW/dWp
// W  = [rho, rho * u, e]
// Wp = [rho, u, p]
void dWpdW(std::vector <double> &M,
           std::vector <double> W,
           int k)
{
    double rho, u;
    rho = W[0 * nx + k];
    u = W[1 * nx + k] / W[0 * nx + k];
    M[0] = 1;
    M[1] = 0;
    M[2] = 0;
    M[3] = -u / rho;
    M[4] = 1 / rho;
    M[5] = 0;
    M[6] = u * u * (gam - 1) / 2;
    M[7] = u * (1 - gam);
    M[8] = gam - 1;
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

// get Q
void WtoQ(std::vector <double> W,
          std::vector <double> &Q,
          std::vector <double> S)
{
    double rho, u, e, p;
    for(int i = 0; i < nx; i++)
    {
        rho = W[0 * nx + i];
        u = W[1 * nx + i] / rho;
        e = W[2 * nx + i];
        p = (gam - 1) * ( e - rho * u * u / 2 );

        Q[1 * nx + i] = p * (S[i + 1] - S[i]);
    }
}
