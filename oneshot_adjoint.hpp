#ifndef ONESHOT_ADJOINT_H
#define ONESHOT_ADJOINT_H

#include<vector>
#include"structures.hpp"
void oneshot_adjoint(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &dx,
	const struct Flow_options &flo_opts,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &initial_design);
#endif
