#ifndef optimizer_h
#define optimizer_h

#include<vector>
#include"structures.h"
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
	const struct Optimization_options &opt_opts,
    VectorXd* const searchD,
    struct Flow_data* const flow_data,
    struct Design* const current_design);

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
	const struct Optimization_options &opt_opts,
	const struct Design &initial_design);
#endif
