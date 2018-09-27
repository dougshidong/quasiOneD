#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include<vector>
#include"structures.h"
void inletBC(
    const Flow_options &flo_opts,
    const double dt0,
	const double dx0,
    struct Flow_data* const flow_data);
void outletBC(
    const Flow_options &flo_opts,
    const double dt0,
	const double dx0,
    struct Flow_data* const flow_data);
#endif
