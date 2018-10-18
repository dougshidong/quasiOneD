#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include<vector>
#include"structures.h"

template<typename dreal>
void inletBC(
    const Flow_options<dreal> &flo_opts,
    const dreal dt0,
	const dreal dx0,
    struct Flow_data<dreal>* const flow_data);

template<typename dreal>
void outletBC(
    const Flow_options<dreal> &flo_opts,
    const dreal dt0,
	const dreal dx0,
    struct Flow_data<dreal>* const flow_data);
#endif
