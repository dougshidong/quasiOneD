#ifndef RESIDUAL_DERIVATIVES_H
#define RESIDUAL_DERIVATIVES_H

#include "structures.h"
#include<vector>

using namespace Eigen;
SparseMatrix<double> dRdW_adolc(
	const struct Flow_options<double> &flo_opts,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    const struct Flow_data<double> &flow_data);

#endif

