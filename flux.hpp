#ifndef FLUX_H
#define FLUX_H

#include "structures.hpp"
#include <vector>

template<typename dreal>
void getFlux(
	const struct Flow_options &flow_options,
    const std::vector<dreal> &W,
    std::vector<dreal> *const fluxes);
#endif
