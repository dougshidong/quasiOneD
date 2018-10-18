#ifndef SECOND_ORDER_FLOW_H
#define SECOND_ORDER_FLOW_H

#include<vector>
#include"structures.h"

void second_order_flow(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &area,
	const struct Flow_options<double> &flo_opts);

#endif
