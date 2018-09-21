#ifndef QUASIONED_H
#define QUASIONED_H

#include<vector>
#include"structures.h"

double isenP(double pt, double M);
double isenT(double Tt, double M);
double quasiOneD(
    const std::vector<double> &x,
    const std::vector<double> &area,
	const Flow_options &flow_options,
    struct Flow_data* const flow_data);
void inletBC(
    const Flow_options &flow_options,
    const double dt0,
	const double dx0,
    struct Flow_data* const flow_data);
void outletBC(
    const Flow_options &flow_options,
    const double dt0,
	const double dx0,
    struct Flow_data* const flow_data);

#endif
