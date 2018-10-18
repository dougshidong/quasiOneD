#ifndef BOUNDARY_GRADIENT_H
#define BOUNDARY_GRADIENT_H
#include <vector>
#include "structures.h"
void dRdW_BC_inlet(
	const struct Flow_options<double> &flow_options,
	const std::vector<double> &W,
    std::vector<double> &dBidWi,
    std::vector<double> &dBidWd);

void dRdW_BC_outlet(
	const struct Flow_options<double> &flow_options,
	const std::vector<double> &W,
    std::vector<double> &dBodWd,
    std::vector<double> &dBodWo);
#endif
