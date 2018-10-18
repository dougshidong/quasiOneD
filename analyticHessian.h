#ifndef ANALYTICHESSIAN_H
#define ANALYTICHESSIAN_H

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
    const struct Flow_options<double> &flo_opts,
    const struct Flow_data<double> &flow_data,
    const struct Optimization_options<double> &opt_opts,
    const struct Design<double> &design);
#endif


