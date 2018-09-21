#ifndef SPLINE_H
#define SPLINE_H
#include "structures.h"
#include <vector>
#include <Eigen/Core>

std::vector<double> evalSpline(
	const int n_control_pts,
	const int spline_degree,
	const std::vector<double> &control_points,
    const std::vector<double> &x,
    const std::vector<double> &dx);

std::vector<double> fit_bspline(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const int n_control_pts,
	const int spline_degree);

Eigen::MatrixXd evalSplineDerivative(
	const int n_control_pts,
	const int spline_degree,
	const std::vector<double> &control_points,
    std::vector<double> x,
    std::vector<double> dx);
#endif
