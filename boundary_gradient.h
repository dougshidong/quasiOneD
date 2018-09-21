#ifndef boundary_gradient_h
#define boundary_gradient_h
#include <vector>
#include "structures.h"
void dRdW_BC_inlet(
	const struct Flow_options &flow_options,
	const struct Flow_data &flow_data,
    std::vector<double> &dBidWi,
    std::vector<double> &dBidWd);

void dRdW_BC_outlet(
	const struct Flow_options &flow_options,
	const struct Flow_data &flow_data,
    std::vector<double> &dBodWd,
    std::vector<double> &dBodWo);
#endif
