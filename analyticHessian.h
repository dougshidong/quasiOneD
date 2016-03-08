#ifndef analyticHessian_h
#define analyticHessian_h

#include<Eigen/Core>

Eigen::MatrixXd getAnalyticHessian(
    std::vector <double> W,
    std::vector <double> S,
    int method);
#endif


