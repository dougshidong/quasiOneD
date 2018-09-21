#include "structures.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "spline.h"

std::vector<double> sinParam(
    const double h, const double t1, const double t2,
    std::vector<double> x,
    std::vector<double> dx)
{
	const double PI = atan(1.0) * 4.0;
	int n_elem = x.size();
    std::vector<double> area(n_elem + 1);
    // Define Area
    for (int i = 1; i < n_elem; i++) {
        double xh = x[i] - dx[i] / 2.0;
        area[i] = 1 - h * pow(sin(PI * pow(xh, t1)), t2);
    }

    area[0] = 1;
    area[n_elem] = 1;

    return area;
}


// Evaluate X
std::vector<double> uniform_x(double a, double b, int n_elem) {
    std::vector<double> x(n_elem);
    double dxConst = (b - a)/n_elem;

    for (int i = 0; i < n_elem; i++)
        x[i] = dxConst/2 + dxConst * i;

    return x;
}

//  Evaluate dx

std::vector<double> eval_dx(std::vector<double> x) {
	int n_elem = x.size();
    std::vector<double> dx(n_elem);

    dx[0] = x[1] - x[0];
    for (int i = 1; i < n_elem - 1; i++) {
        dx[i] =  (x[i] - x[i - 1])/2  +  (x[i + 1] - x[i])/2 ;
    }
    dx[n_elem - 1] = x[x.size() - 1] - x[x.size() - 2];

    return dx;
}

std::vector<double> evalS(
	const struct Design &design,
    const std::vector<double> &x,
    const std::vector<double> &dx)
{
	int n_elem = x.size();
    std::vector<double> area(n_elem + 1);
    if (design.parametrization == 0) {
		abort();
    }
    else if (design.parametrization == 1) {
        area = sinParam(design.h, design.t1, design.t2, x, dx);\
    }
    else if (design.parametrization == 2) {
		int n_control_pts = design.n_design_variables + 2;
		int spline_degree = design.spline_degree;
        std::vector<double> control_points;
		control_points.push_back(1); // Clamped
		control_points.insert(control_points.end(), design.design_variables.begin(), design.design_variables.end());
		control_points.push_back(1); // Clamped

        area = evalSpline(n_control_pts, spline_degree, control_points, x, dx);
    }
    return area;
}

