#ifndef FLUX_H
#define FLUX_H

#include "structures.h"
#include <vector>

void getFlux(
	const struct Flow_options &flow_options,
    const std::vector<double> &W,
    std::vector<double> &fluxes);
#endif
