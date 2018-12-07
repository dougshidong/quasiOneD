#ifndef QUASIONED_H
#define QUASIONED_H

#include<vector>
#include"structures.hpp"

template<typename dreal>
int quasiOneD(
	const bool restart,
    const std::vector<dreal> &x,
    const std::vector<dreal> &area,
	const Flow_options &flow_options,
    class Flow_data<dreal>* const flow_data);
#endif
