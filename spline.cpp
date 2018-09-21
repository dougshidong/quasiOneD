#include "structures.h"
#include <iostream>
#include <math.h>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SVD>

using namespace Eigen;

std::vector<double> getKnots(
	const int nknots,
	const int spline_degree);
double getbij(double x, int i, int j, std::vector<double> t);

std::vector<double> evalSpline(
	const int n_control_pts,
	const int spline_degree,
	const std::vector<double> &control_points,
    const std::vector<double> &x,
    const std::vector<double> &dx)
{
	int n_elem = x.size();
    // Get Knot Vector
    int nknots = n_control_pts + spline_degree + 1;
    std::vector<double> knots(nknots);
    knots = getKnots(nknots, spline_degree);

    // Define Area at i+1/2
    std::vector<double> area(n_elem + 1, 0);
    for (int Si = 1; Si < n_elem; Si++)
    {
        double xh = x[Si] - dx[Si] / 2.0;
        for (int ictl = 0; ictl < n_control_pts; ictl++)
        {
            area[Si] += control_points[ictl] * getbij(xh, ictl, spline_degree, knots);
        }
    }
    area[0] = 1;
    area[n_elem] = 1;

    return area;
}


MatrixXd evalSplineDerivative( // now returns dSdCtl instead of dSdDesign (nctl = ndes+2)
	const int n_control_pts,
	const int spline_degree,
	const std::vector<double> &control_points,
    std::vector<double> x,
    std::vector<double> dx)
{
	int n_elem = x.size();
    MatrixXd dSdCtl(n_elem + 1, n_control_pts);
    dSdCtl.setZero();

    // Get Knot Vector
    int nknots = n_control_pts + spline_degree + 1;
    std::vector<double> knots(nknots);
    knots = getKnots(nknots, spline_degree);

    // Define Area at i+1/2
	double xh;
    for (int Si = 0; Si < n_elem + 1; Si++) {
        if (Si < n_elem) xh = x[Si] - dx[Si] / 2.0;
        else xh = 1.0;
        for (int ictl = 0; ictl < n_control_pts; ictl++) {// Not including the inlet/outlet
            dSdCtl(Si, ictl) += getbij(xh, ictl, spline_degree, knots);
        }
    }
    return dSdCtl;
}
MatrixXd eval_dArea_dDesign(
	const struct Design &design,
    std::vector<double> x,
    std::vector<double> dx)
{
	int n_elem = x.size();

	std::vector<double> control_points;
	control_points.push_back(1); // Clamped
	control_points.insert(control_points.end(), design.design_variables.begin(), design.design_variables.end());
	control_points.push_back(1); // Clamped
	int n_control_pts = design.n_design_variables+2;
	int spline_degree = design.spline_degree;

    // Get Knot Vector
    int nknots = n_control_pts + spline_degree + 1;
    std::vector<double> knots(nknots);
    knots = getKnots(nknots, spline_degree);

	MatrixXd dArea_dSpline = evalSplineDerivative(n_control_pts, spline_degree, control_points, x, dx);

	int block_start_i = 0;
	int block_start_j = 0;
	int i_size = n_elem + 1;
	int j_size = design.n_design_variables;
    return dArea_dSpline.block(block_start_i, block_start_j, i_size, j_size); // Do not return the endpoints
}

std::vector<double> fit_bspline(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const int n_control_pts,
	const int spline_degree)
{
	int n_elem = x.size();
    // Get Knot Vector
    int nknots = n_control_pts + spline_degree + 1;
    std::vector<double> knots(nknots);
    knots = getKnots(nknots, spline_degree);

    VectorXd ctl_eig(n_control_pts), s_eig(n_elem + 1);
    MatrixXd A(n_elem + 1, n_control_pts);
    A.setZero();
    
	double xend = 1.0;
	double xh;
    for (int iface = 0; iface < n_elem + 1; iface++) {
        if (iface < n_elem) xh = x[iface] - dx[iface] / 2.0;
        else xh = xend;
        s_eig(iface) = area[iface];
        for (int ictl = 0; ictl < n_control_pts; ictl++)
        {
            A(iface, ictl) = getbij(xh, ictl, spline_degree, knots);
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

double getbij(double x, int i, int j, std::vector<double> t) {
    if (j==0) {
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

std::vector<double> getKnots(
	const int nknots,
	const int spline_degree)
{
	double grid_xstart = 0.0;
	double grid_xend   = 1.0;
    std::vector<double> knots(nknots);
    int nb_outer = 2 * (spline_degree + 1);
    int nb_inner = nknots - nb_outer;
    double eps = 2e-15; // Allow Spline Definition at End Point
    // Clamped Open-Ended
    for (int iknot = 0; iknot < spline_degree + 1; iknot++)
    {
        knots[iknot] = grid_xstart;
        knots[nknots - iknot - 1] = grid_xend + eps;
    }
    // Uniform Knot Vector
    double knot_dx = (grid_xend + eps - grid_xstart) / (nb_inner + 1);
    for (int iknot = 1; iknot < nb_inner+1; iknot++)
    {
        knots[iknot + spline_degree] = iknot * knot_dx;
    }
    return knots;
}
