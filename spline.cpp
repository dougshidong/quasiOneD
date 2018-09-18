#include <iostream>
#include <math.h>
#include <vector>
#include "globals.h"
#include <Eigen/Core>
#include <Eigen/SVD>

using namespace Eigen;

std::vector<double> getKnots(int nknots);
double getbij(double x, int i, int j, std::vector<double> t);

std::vector<double> evalSpline(
    std::vector<double> ctlpts,
    std::vector<double> x,
    std::vector<double> dx)
{
    // Get Knot Vector
    int nknots = n_control_pts + spline_degree + 1;
    std::vector<double> knots(nknots);
    knots = getKnots(nknots);

    // Define Area at i+1/2
    std::vector<double> area(n_elem + 1, 0);
    double xh;
    for (int Si = 1; Si < n_elem; Si++)
    {
        xh = x[Si] - dx[Si] / 2.0;
        for (int ictl = 0; ictl < n_control_pts; ictl++)
        {
            area[Si] += ctlpts[ictl] * getbij(xh, ictl, spline_degree, knots);
        }
    }
    area[0] = 1;
    area[n_elem] = 1;

    return area;
}

MatrixXd evalSplineDerivative(
    std::vector<double> x,
    std::vector<double> dx)
{
    MatrixXd dSdCtl(n_elem + 1, nDesVar);
    dSdCtl.setZero();

    // Get Knot Vector
    int nknots = n_control_pts + spline_degree + 1;
    std::vector<double> knots(nknots);
    knots = getKnots(nknots);

    // Define Area at i+1/2
    double xh;
    for (int Si = 0; Si < n_elem + 1; Si++)
    {
        xh = x[Si] - dx[Si] / 2.0;
        for (int ictl = 1; ictl < n_control_pts - 1; ictl++) // Not including the inlet/outlet
        {
            dSdCtl(Si, ictl - 1) += getbij(xh, ictl, spline_degree, knots);
        }
    }
    return dSdCtl;
}

std::vector<double> getCtlpts(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> area)
{
    // Get Knot Vector
    int nknots = n_control_pts + spline_degree + 1;
    std::vector<double> knots(nknots);
    knots = getKnots(nknots);

    VectorXd ctl_eig(n_control_pts), s_eig(n_elem + 1);
    MatrixXd A(n_elem + 1, n_control_pts);
    A.setZero();
    
    double xh;
    for (int Si = 0; Si < n_elem + 1; Si++)
    {
        if (Si < n_elem) xh = x[Si] - dx[Si] / 2.0;
        else xh = b_geom;
        s_eig(Si) = area[Si];
        for (int ictl = 0; ictl < n_control_pts; ictl++)
        {
            A(Si, ictl) = getbij(xh, ictl, spline_degree, knots);
        }
    }
    ctl_eig = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(s_eig);
    
    std::vector<double> ctlpts(n_control_pts);
    for (int ictl = 0; ictl < n_control_pts; ictl++)
    {
        ctlpts[ictl] = ctl_eig(ictl);
    }
    ctlpts[0] = area[0];
    ctlpts[n_control_pts - 1] = area[n_elem];

    return ctlpts;
}

double getbij(double x, int i, int j, std::vector<double> t)
{
    if (j==0)
    {
        if (t[i] <= x && x < t[i+1]) return 1;
        else return 0;
    }

    double h = getbij(x, i,   j-1, t);
    double k = getbij(x, i+1, j-1, t);

    double bij = 0;

    if (h!=0) bij += (x        - t[i]) / (t[i+j]   - t[i]  ) * h;
    if (k!=0) bij += (t[i+j+1] - x   ) / (t[i+j+1] - t[i+1]) * k;

    return bij;
}

std::vector<double> getKnots(int nknots)
{
    std::vector<double> knots(nknots);
    int nb_outer = 2 * (spline_degree + 1);
    int nb_inner = nknots - nb_outer;
    double eps = 2e-15; // Allow Spline Definition at End Point
    // Clamped Open-Ended
    for (int iknot = 0; iknot < spline_degree + 1; iknot++)
    {
        knots[iknot] = a_geom;
        knots[nknots - iknot - 1] = b_geom + eps;
    }
    // Uniform Knot Vector
    double knot_dx = (b_geom + eps - a_geom) / (nb_inner + 1);
    for (int iknot = 1; iknot < nb_inner+1; iknot++)
    {
        knots[iknot + spline_degree] = iknot * knot_dx;
    }
    return knots;
}
