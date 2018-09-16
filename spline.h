#ifndef spline_h
#define spline_h
#include <vector>
#include <Eigen/Core>
std::vector<double> evalSpline(
    std::vector<double> geom,
    std::vector<double> x,
    std::vector<double> dx);

std::vector<double> getCtlpts(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> area);

Eigen::MatrixXd evalSplineDerivative(
    std::vector<double> x,
    std::vector<double> dx);
#endif
