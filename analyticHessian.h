#ifndef analyticHessian_h
#define analyticHessian_h

#include<Eigen/Core>
#include"structures.h"
// case 0: directDirectHessian
// case 1: adjointDirectHessian
// case 2: adjointAdjointHessian
// case 3: directAdjointHessian
Eigen::MatrixXd getAnalyticHessian(
    const int hessian_type,
    const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
    const struct Flow_options &flo_opts,
    const struct Flow_data &flow_data,
    const struct Optimization_options &opt_opts,
    const struct Design &design);
#endif


