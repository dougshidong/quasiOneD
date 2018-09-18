#ifndef analyticHessian_h
#define analyticHessian_h

#include<Eigen/Core>

Eigen::MatrixXd getAnalyticHessian(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar,
    int method);
#endif


