// Convert W to F or primitive variables

#include<vector>
#include<math.h>
#include"globals.h"


// Get primitive variables
void WtoP(std::vector <std::complex<double> > W,
          std::vector <std::complex<double> > &rho,
          std::vector <std::complex<double> > &u,
          std::vector <std::complex<double> > &e)
{
    for(int i = 0; i < nx; i++)
    {
        rho[i] = W[i * 3 + 0];
        u[i] = W[i * 3 + 1] / rho[i];
        e[i] = W[i * 3 + 2];
    }
}

// Get more primitive variables
void WtoP(std::vector <std::complex<double> > W,
          std::vector <std::complex<double> > &rho,
          std::vector <std::complex<double> > &u,
          std::vector <std::complex<double> > &e,
          std::vector <std::complex<double> > &p,
          std::vector <std::complex<double> > &c,
          std::vector <std::complex<double> > &T)
{
    for(int i = 0; i < nx; i++)
    {
        rho[i] = W[i * 3 + 0];
        u[i] = W[i * 3 + 1] / rho[i];
        e[i] = W[i * 3 + 2];
        p[i] = (gam - 1) * ( e[i] - rho[i] * u[i] * u[i] / 2.0 );
        c[i] = sqrt( gam * p[i] / rho[i] );
        T[i] = p[i] / (rho[i] * R);
    }
}
// Inverse[dW/dWp] or dWp/dW
// W  = [rho, rho * u, e]
// Wp = [rho, u, p]
void dWpdW(std::vector <std::complex<double> > &M,
           std::vector <std::complex<double> > W,
           int i)
{
    std::complex<double> rho, u;
    rho = W[i * 3 + 0];
    u = W[i * 3 + 1] / rho;
    M[0] = 1.0;
    M[1] = 0.0;
    M[2] = 0.0;
    M[3] = -u / rho;
    M[4] = 1.0 / rho;
    M[5] = 0.0;
    M[6] = u * u * (gam - 1.0) / 2.0;
    M[7] = u * (1.0 - gam);
    M[8] = gam - 1.0;
}
// dW/dWp
// W  = [rho, rho * u, e]
// Wp = [rho, u, p]
void dWdWp(std::vector <std::complex<double> > &M,
           std::vector <std::complex<double> > W,
           int i)
{
    std::complex<double> rho, u;
    rho = W[i * 3 + 0];
    u = W[i * 3 + 1] / rho;
    M[0] = 1.0;
    M[1] = 0.0;
    M[2] = 0.0;
    M[3] = u;
    M[4] = rho;
    M[5] = 0.0;
    M[6] = u * u / 2.0;
    M[7] = u * rho;
    M[8] = 1.0 / (gam - 1.0);
}
// Get F
void WtoF(std::vector <std::complex<double> > W,
          std::vector <std::complex<double> > &F)
{
    std::complex<double> w1, w2, w3;
    for(int i = 0; i < nx; i++)
    {
        w1 = W[i * 3 + 0];
        w2 = W[i * 3 + 1];
        w3 = W[i * 3 + 2];
        F[i * 3 + 0] = w2;
        F[i * 3 + 1] = w2 * w2 / w1 + (gam - 1) * ( w3 - w2 * w2 / (2.0 * w1) );
        F[i * 3 + 2] = ( w3 + (gam - 1) * (w3 - w2 * w2 / (2.0 * w1)) ) * w2 / w1;
    }
}

// get Q
void WtoQ(std::vector <std::complex<double> > W,
          std::vector <std::complex<double> > &Q,
          std::vector <std::complex<double> > S)
{
    std::complex<double> rho, u, e, p;
    for(int i = 0; i < nx; i++)
    {
        rho = W[i * 3 + 0];
        u = W[i * 3 + 1] / rho;
        e = W[i * 3 + 2];
        p = (gam - 1) * ( e - rho * u * u / 2.0 );

        Q[i * 3 + 0] = 0;
        Q[i * 3 + 1] = p * (S[i + 1] - S[i]);
        Q[i * 3 + 2] = 0;
    }
}
