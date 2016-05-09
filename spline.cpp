#include <iostream>
#include <math.h>
#include <vector>
#include "globals.h"
#include <Eigen/Core>
#include <Eigen/SVD>

using namespace Eigen;

std::vector <double> getKnots(int nknots);
double getbij(double x, int i, int j, std::vector <double> t);

std::vector <double> evalSpline(
    std::vector <double> ctlpts,
    std::vector <double> x,
    std::vector <double> dx)
{
    // Get Knot Vector
    int nknots = nctl + spline_degree + 1;
    std::vector <double> knots(nknots);
    knots = getKnots(nknots);

    // Define Area at i+1/2
    std::vector <double> S(nx + 1, 0);
    double xh;
    for(int Si = 1; Si < nx; Si++)
    {
        xh = x[Si] - dx[Si] / 2.0;
        for(int ictl = 0; ictl < nctl; ictl++)
        {
            S[Si] += ctlpts[ictl] * getbij(xh, ictl, spline_degree, knots);
        }
    }
    S[0] = 1;
    S[nx] = 1;

    return S;
}

MatrixXd evalSplineDerivative(
    std::vector <double> x,
    std::vector <double> dx)
{
    MatrixXd dSdCtl(nx + 1, nDesVar);
    dSdCtl.setZero();

    // Get Knot Vector
    int nknots = nctl + spline_degree + 1;
    std::vector <double> knots(nknots);
    knots = getKnots(nknots);

    // Define Area at i+1/2
    double xh;
    for(int Si = 0; Si < nx + 1; Si++)
    {
        xh = x[Si] - dx[Si] / 2.0;
        for(int ictl = 1; ictl < nctl - 1; ictl++) // Not including the inlet/outlet
        {
            dSdCtl(Si, ictl - 1) += getbij(xh, ictl, spline_degree, knots);
        }
    }
    return dSdCtl;
}

std::vector <double> getCtlpts(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S)
{
    // Get Knot Vector
    int nknots = nctl + spline_degree + 1;
    std::vector <double> knots(nknots);
    knots = getKnots(nknots);

    VectorXd ctl_eig(nctl), s_eig(nx + 1);
    MatrixXd A(nx + 1, nctl);
    A.setZero();
    
    double xh;
    for(int Si = 0; Si < nx + 1; Si++)
    {
        if(Si < nx) xh = x[Si] - dx[Si] / 2.0;
        else xh = b_geom;
        s_eig(Si) = S[Si];
        for(int ictl = 0; ictl < nctl; ictl++)
        {
            A(Si, ictl) = getbij(xh, ictl, spline_degree, knots);
        }
    }
    ctl_eig = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(s_eig);
    
    std::vector <double> ctlpts(nctl);
    for(int ictl = 0; ictl < nctl; ictl++)
    {
        ctlpts[ictl] = ctl_eig(ictl);
    }
    ctlpts[0] = S[0];
    ctlpts[nctl - 1] = S[nx];

    return ctlpts;
}

double getbij(double x, int i, int j, std::vector <double> t)
{
    if(j==0)
    {
        if(t[i] <= x && x < t[i+1]) return 1;
        else return 0;
    }

    double h = getbij(x, i,   j-1, t);
    double k = getbij(x, i+1, j-1, t);

    double bij = 0;

    if(h!=0) bij += (x        - t[i]) / (t[i+j]   - t[i]  ) * h;
    if(k!=0) bij += (t[i+j+1] - x   ) / (t[i+j+1] - t[i+1]) * k;

    return bij;
}

std::vector <double> getKnots(int nknots)
{
    std::vector <double> knots(nknots);
    int nb_outer = 2 * (spline_degree + 1);
    int nb_inner = nknots - nb_outer;
    double eps = 2e-15; // Allow Spline Definition at End Point
    // Clamped Open-Ended
    for(int iknot = 0; iknot < spline_degree + 1; iknot++)
    {
        knots[iknot] = a_geom;
        knots[nknots - iknot - 1] = b_geom + eps;
    }
    // Uniform Knot Vector
    double knot_dx = (b_geom + eps - a_geom) / (nb_inner + 1);
    for(int iknot = 1; iknot < nb_inner+1; iknot++)
    {
        knots[iknot + spline_degree] = iknot * knot_dx;
    }
    return knots;
}
