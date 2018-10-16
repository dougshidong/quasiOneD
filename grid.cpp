#include "structures.h"
#include <iostream>
#include <math.h>
#include <vector>
#include "spline.h"

std::vector<double> sinParam(
    const double h, const double t1, const double t2,
    std::vector<double> x)
{
	const double PI = atan(1.0) * 4.0;
	const int n_face = x.size();
    std::vector<double> area(n_face);
    // Define Area
    for (int i = 0; i < n_face; i++) {
        area[i] = 1 - h * pow(sin(PI * pow(x[i], t1)), t2);
    }

    return area;
}


// Evaluate X
std::vector<double> uniform_x(double a, double b, int n_elem) {
    const int n_face = n_elem+1;
    std::vector<double> x(n_face);
    const double dxConst = (b - a)/n_elem;
    for (int i = 0; i < n_face; i++) {
        x[i] = dxConst * i;
        //x[i] = dxConst/2 + dxConst * i;
    }
    return x;
}

//  Evaluate dx

std::vector<double> eval_dx(std::vector<double> x) {
	const int n_face = x.size();
	const int n_dx = n_face+1;
    std::vector<double> dx(n_dx);
    for (int i = 0; i < n_face-1; i++) {
        dx[i+1] =  x[i+1] - x[i];
    }
    dx[0] = dx[1];
    dx[n_face] = dx[n_face-1];
    return dx;
}

std::vector<double> evalS(
	const struct Design &design,
    const std::vector<double> &x,
    const std::vector<double> &dx)
{
	int n_face = x.size();
    std::vector<double> area(n_face);
    if (design.parametrization == 0) {
		abort();
    }
    else if (design.parametrization == 1) {
        area = sinParam(design.h, design.t1, design.t2, x);\
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

