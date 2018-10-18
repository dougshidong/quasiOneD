#ifndef FLUX_H
#define FLUX_H

#include "structures.h"
#include <vector>

template<typename dreal>
void getFlux(
	const struct Flow_options<dreal> &flow_options,
    const std::vector<dreal> &W,
    std::vector<dreal> *const fluxes);
#endif
