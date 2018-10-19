#ifndef SPLINE_H
#define SPLINE_H
#include "structures.hpp"
#include <vector>
#include <Eigen/Core>

template<typename dreal>
std::vector<dreal> evalSpline(
	const int n_control_pts,
	const int spline_degree,
	const std::vector<dreal> &control_points,
    const std::vector<dreal> &x,
    const std::vector<dreal> &dx);

template<typename dreal>
std::vector<dreal> fit_bspline(
    const std::vector<dreal> &x,
    const std::vector<dreal> &dx,
    const std::vector<dreal> &area,
	const int n_control_pts,
	const int spline_degree);

Eigen::MatrixXd evalSplineDerivative(
	const int n_control_pts,
	const int spline_degree,
	const std::vector<double> &control_points,
    std::vector<double> x,
    std::vector<double> dx);
#endif
