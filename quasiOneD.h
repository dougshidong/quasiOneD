#ifndef QUASIONED_H
#define QUASIONED_H

#include<vector>
#include"structures.h"

template<typename dreal>
dreal isenP(const dreal gam, const dreal pt, const dreal M);

template<typename dreal>
dreal isenT(const dreal gam, const dreal Tt, const dreal M);

template<typename dreal>
int quasiOneD(
    const std::vector<dreal> &x,
    const std::vector<dreal> &area,
	const Flow_options<dreal> &flow_options,
    struct Flow_data<dreal>* const flow_data);
#endif
