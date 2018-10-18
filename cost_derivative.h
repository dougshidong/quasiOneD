#ifndef COST_DERIVATIVE_H
#define COST_DERIVATIVE_H
#include"structures.h"
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<vector>

using namespace Eigen;
VectorXd evaldCostdW(
	const struct Optimization_options<double> &opt_opts,
	const struct Flow_options<double> &flo_opts,
	const std::vector<double> &W,
	const std::vector<double> &dx);
SparseMatrix<double> evaldCostdWdW(
	const struct Optimization_options<double> &opt_opts,
	const struct Flow_options<double> &flo_opts,
	std::vector<double> W,
	std::vector<double> dx);

VectorXd evaldCostdArea(const int n_elem);
MatrixXd evalddCostdAreadArea(const int n_elem);
MatrixXd evalddCostdWdArea(const int n_elem);
#endif
