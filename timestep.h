#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "structures.h"
#include<vector>

void initializeTimeStep(int n_elem);
void getDomainResi( 
	const struct Flow_options &flo_opts,
	const std::vector<double> &area,
	const std::vector<double> &W,
	std::vector<double> *const fluxes,
	std::vector<double> *const residual);
void stepInTime(
	const struct Flow_options &flow_options,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    struct Flow_data* const flow_data);

#endif
