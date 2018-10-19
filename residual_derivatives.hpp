#ifndef RESIDUAL_DERIVATIVES_H
#define RESIDUAL_DERIVATIVES_H

#include "structures.hpp"
#include<vector>

Eigen::SparseMatrix<double> eval_dRdW_dRdX_adolc(
	const struct Flow_options &flo_opts,
    const std::vector<double> &area,
    const class Flow_data<double> &flow_data);

#endif

