#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "structures.h"
#include<vector>

template<typename dreal>
void getDomainResi( 
	const struct Flow_options<dreal> &flo_opts,
	const std::vector<dreal> &area,
	const std::vector<dreal> &W,
	std::vector<dreal>* const fluxes,
	std::vector<dreal>* const residual);

template<typename dreal>
void stepInTime(
	const struct Flow_options<dreal> &flow_options,
    const std::vector<dreal> &area,
    const std::vector<dreal> &dx,
    struct Flow_data<dreal>* const flow_data);

#endif
