#ifndef HESSIANINLET_H
#define HESSIANINLET_H

#include<vector>
#include<Eigen/Core>
#include"structures.h"

using namespace Eigen;

void HessianInlet(
    const struct Flow_options<double> &flo_opts,
    const std::vector<double> &W,
    std::vector <MatrixXd> &ddRindWdW);
#endif
