#ifndef optimizer_h
#define optimizer_h

#include<vector>
#include"structures.h"
void optimizer(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &dx,
	const struct Flow_options &flo_opts,
	const struct Optimization_options &opt_opts,
	const struct Design &initial_design);
#endif
