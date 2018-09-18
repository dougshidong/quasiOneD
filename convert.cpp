// Convert W to F or primitive variables

#include<vector>
#include<Eigen/Core>
#include<math.h>
#include"globals.h"


// Get pressure
void getp(
    const std::vector<double> &W,
    std::vector<double> &p)
{
    for (int i = 0; i < n_elem; i++)
    {
        p[i] = (gam - 1.0) * ( W[i * 3 + 2] - (pow(W[i * 3 + 1], 2.0) / W[i * 3 + 0]) / 2.0 );
    }
}

// Get primitive variables
void WtoP(
    const std::vector<double> &W,
    std::vector<double> &rho,
    std::vector<double> &u,
    std::vector<double> &e)
{
    for (int i = 0; i < n_elem; i++)
    {
        rho[i] = W[i * 3 + 0];
        u[i] = W[i * 3 + 1] / rho[i];
        e[i] = W[i * 3 + 2];
    }
}

// Get other primitive variables
void WtoP2(
    const std::vector<double> &W,
    std::vector<double> &rho,
    std::vector<double> &u,
    std::vector<double> &p)
{
    for (int i = 0; i < n_elem; i++)
    {
        rho[i] = W[i * 3 + 0];
        u[i] = W[i * 3 + 1] / rho[i];
        p[i] = (gam - 1.0) * ( W[i * 3 + 2] - (pow(W[i * 3 + 1], 2.0) / W[i * 3 + 0]) / 2.0 );
    }
}

// Given rho, u p, get W
void PtoW(
    std::vector<double> &W,
    const std::vector<double> &Wp)
{
    for (int i = 0; i < n_elem; i++)
    {
        W[i * 3 + 0] = Wp[i * 3 + 0];
        W[i * 3 + 1] = Wp[i * 3 + 0] * Wp[i * 3 + 1];
        W[i * 3 + 2] = Wp[i * 3 + 2] / (gam - 1.0) + 
                       (pow(Wp[i * 3 + 1], 2.0) * Wp[i * 3 + 0]) / 2.0;
    }
}

// Get more primitive variables
void WtoP(
    const std::vector<double> &W,
    std::vector<double> &rho,
    std::vector<double> &u,
    std::vector<double> &e,
    std::vector<double> &p,
    std::vector<double> &c,
    std::vector<double> &T)
{
    for (int i = 0; i < n_elem; i++)
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
void dWpdW(
    std::vector<double> &M,
    const std::vector<double> &W,
    int i)
{
    double rho, u;
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

Eigen::MatrixXd dWpdW(
    const std::vector<double> &W,
    int i)
{
    Eigen::MatrixXd M(3, 3);
    double rho, u;
    rho = W[i * 3 + 0];
    u = W[i * 3 + 1] / rho;
    M(0,0) = 1.0;
    M(0,1) = 0.0;
    M(0,2) = 0.0;
    M(1,0) = -u / rho;
    M(1,1) = 1.0 / rho;
    M(1,2) = 0.0;
    M(2,0) = u * u * (gam - 1.0) / 2.0;
    M(2,1) = u * (1.0 - gam);
    M(2,2) = gam - 1.0;
    return M;
}



// dW/dWp
// W  = [rho, rho * u, e]
// Wp = [rho, u, p]
void dWdWp(
    std::vector<double> &M,
    const std::vector<double> &W,
    int i)
{
    double rho, u;
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
void WtoF(
    const std::vector<double> &W,
    std::vector<double> &F)
{
    double w1, w2, w3;
    for (int i = 0; i < n_elem; i++)
    {
        w1 = W[i * 3 + 0];
        w2 = W[i * 3 + 1];
        w3 = W[i * 3 + 2];
        F[i * 3 + 0] = w2;
        F[i * 3 + 1] = w2 * w2 / w1 + (gam - 1) * ( w3 - w2 * w2 / (2 * w1) );
        F[i * 3 + 2] = ( w3 + (gam - 1) * (w3 - w2 * w2 / (2 * w1)) ) * w2 / w1;
    }
}

// get Q
void WtoQ(
    const std::vector<double> &W,
    std::vector<double> &Q,
    const std::vector<double> &area)
{
    std::vector<double> p(n_elem);
    getp(W, p);
    for (int i = 0; i < n_elem; i++)
    {
        Q[i * 3 + 0] = 0;
        Q[i * 3 + 1] = p[i] * (area[i + 1] - area[i]);
        Q[i * 3 + 2] = 0;
    }
}
