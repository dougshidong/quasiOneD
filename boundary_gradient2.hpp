#ifndef BOUNDARY_GRADIENT2_H
#define BOUNDARY_GRADIENT2_H

#include<vector>
#include"structures.hpp"

template<typename dreal>
Eigen::Matrix<dreal,3,3> inletBC_gradient(
    const Flow_options &flo_opts,
    class Flow_data<dreal>* const flow_data);

template<typename dreal>
Eigen::Matrix<dreal,3,3> outletBC_gradient(
    const Flow_options &flo_opts,
    class Flow_data<dreal>* const flow_data);
#endif

