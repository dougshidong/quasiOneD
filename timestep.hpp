#ifndef TIMESTEP_H
#define TIMESTEP_H

#include "structures.hpp"
#include<vector>

template<typename dreal>
void getDomainResi( 
	const struct Flow_options &flo_opts,
	const std::vector<dreal> &area,
    class Flow_data<dreal>* const flow_data);

template<typename dreal>
void stepInTime(
	const struct Flow_options &flow_options,
    const std::vector<dreal> &area,
    const std::vector<dreal> &dx,
    class Flow_data<dreal>* const flow_data);

#endif
