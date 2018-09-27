#ifndef cost_derivatives_h
#define cost_derivatives_h
#include<Eigen/Core>
#include<vector>

using namespace Eigen;
VectorXd evaldCostdW(
	const Optimization_options &opt_opts,
	const Flow_options &flo_opts,
	const std::vector<double> &W,
	const std::vector<double> &dx);
SparseMatrix<double> evaldCostdWdW(
	const Optimization_options &opt_opts,
	const Flow_options &flo_opts,
	std::vector<double> W,
	std::vector<double> dx);

VectorXd evaldCostdArea(const int n_elem);
MatrixXd evalddCostdAreadArea(const int n_elem);
MatrixXd evalddCostdWdArea(const int n_elem);
#endif
