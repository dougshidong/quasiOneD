#ifndef FIXED_POINT_DERIVATIVES_H
#define FIXED_POINT_DERIVATIVES_H

#include "structures.hpp"
#include<vector>

void eval_dGdW_dGdX_adolc(
    const std::vector<double> &x,
	const struct Flow_options &flo_opts,
    const class Flow_data<double> &flow_data,
	const struct Design<double> &design,
	SparseMatrix<double> *const dGdW,
	SparseMatrix<double> *const dGdX);

#endif
