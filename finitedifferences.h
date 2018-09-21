#ifndef FINITEDIFFERENCES_H
#define FINITEDIFFERENCES_H

#include "structures.h"
#include <vector>
#include <Eigen/Core>

MatrixXd hessian_central_gradient(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flow_options,
	const struct Flow_data &flow_data,
	const struct Optimization_options &opt_opts,
	const struct Design &design,
    double pert)
MatrixXd hessian_central(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options &flow_options,
	const struct Flow_data &flow_data,
	const struct Optimization_options &opt_opts,
	const struct Design &design,
    double pert)
#endif
