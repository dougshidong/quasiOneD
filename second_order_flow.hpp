#ifndef SECOND_ORDER_FLOW_H
#define SECOND_ORDER_FLOW_H

#include<vector>
#include"structures.hpp"

void second_order_flow(
	const struct Constants &constants,
    const std::vector<double> &x,
	const std::vector<double> &area,
	const struct Flow_options &flo_opts);

#endif
