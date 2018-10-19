#ifndef ONESHOT_DWDX_H
#define ONESHOT_DWDX_H

#include<vector>
#include"structures.hpp"
void oneshot_dwdx(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &dx,
	const struct Flow_options &flo_opts,
	const struct Optimization_options<double> &opt_opts,
	const struct Design<double> &initial_design);
#endif
