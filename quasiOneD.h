#ifndef QUASIONED_H
#define QUASIONED_H

#include<vector>
#include"structures.h"

double isenP(const double gam, const double pt, const double M);
double isenT(const double gam,double Tt, const double M);
double quasiOneD(
    const std::vector<double> &x,
    const std::vector<double> &area,
	const Flow_options &flow_options,
    struct Flow_data* const flow_data);
#endif
