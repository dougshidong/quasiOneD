#ifndef FINITEDIFFERENCES_H
#define FINITEDIFFERENCES_H

#include "structures.h"
#include <vector>
#include <Eigen/Core>

using namespace Eigen;
MatrixXd hessian_central_gradient(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options<double> &flow_options,
	const struct Flow_data<double> &flow_data,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &design,
    double pert);
MatrixXd hessian_central(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options<double> &flow_options,
	const struct Flow_data<double> &flow_data,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &design,
    double pert);
#endif
