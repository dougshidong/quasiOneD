#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include<vector>
#include"structures.hpp"
#include<Eigen/Dense>
using namespace Eigen;
double linesearch_backtrack_unconstrained(
    const double initial_alpha,
    const std::vector<double> x,
    const std::vector<double> dx,
    const VectorXd &pk,
    const VectorXd &gradient,
    const double current_cost,
	const struct Flow_options &flo_opts,
	const struct Optimization_options<double> &opt_opts,
    VectorXd* const searchD,
    class Flow_data<double>* const flow_data,
    struct Design<double>* const current_design);

MatrixXd BFGS(
    const MatrixXd &oldH,
    const VectorXd &oldg,
    const VectorXd &currentg,
    const VectorXd &searchD);

double checkCond(MatrixXd H);
MatrixXd invertHessian(MatrixXd H);
LLT<MatrixXd> checkPosDef(MatrixXd H);

VectorXd implicitSmoothing(VectorXd gradient, double epsilon);
void optimizer(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &dx,
	const struct Flow_options &flo_opts,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &initial_design);
#endif
