#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "structures.h"
#include<vector>

void initializeTimeStep(int n_elem);
void getDomainResi( 
	const struct Flow_options &flow_options,
	const std::vector<double> &area,
	const std::vector<double> &W,
	const std::vector<double> &fluxes,
	std::vector<double> &residual);
void stepInTime(
	const struct Flow_options &flow_options,
    const std::vector<double> &area,
    const std::vector<double> &dx,
    const std::vector<double> &dt,
    struct Flow_data* const flow_data);

#endif
