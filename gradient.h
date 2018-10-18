#ifndef GRADIENT_H
#define GRADIENT_H
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <vector>
#include "structures.h"

using namespace Eigen;

VectorXd getGradient(
	const int gradient_type,
	const int cost_function,
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options<double> &flow_options,
	const struct Flow_data<double> &flow_data,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &design);

MatrixXd evaldWdDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &area,
	const struct Flow_options<double> &flow_options,
	const struct Flow_data<double> &flow_data,
	const struct Design<double> &design);

#endif
