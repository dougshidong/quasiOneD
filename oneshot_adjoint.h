#ifndef oneshot_dwdx_h
#define oneshot_dwdx_h

#include<vector>
#include"structures.h"
void oneshot_dwdx(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &dx,
	const struct Flow_options &flo_opts,
	const struct Optimization_options &opt_opts,
	const struct Design &initial_design);
#endif
