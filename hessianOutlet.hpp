#ifndef HESSIANOUTLET_H
#define HESSIANOUTLET_H

#include<vector>
#include<Eigen/Core>
#include"structures.hpp"

using namespace Eigen;

void HessianOutlet(
    const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector <MatrixXd> &ddRindWdW);
#endif
