#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include<vector>
#include"structures.hpp"

template<typename dreal>
void inletBC(
    const Flow_options &flo_opts,
    const dreal dt0,
	const dreal dx0,
    class Flow_data<dreal>* const flow_data);

template<typename dreal>
void outletBC(
    const Flow_options &flo_opts,
    const dreal dt0,
	const dreal dx0,
    class Flow_data<dreal>* const flow_data);
#endif
