
#ifndef hessianInlet_h
#define hessianInlet_h

#include<vector>
#include<Eigen/Core>
#include"structures.h"

using namespace Eigen;

void HessianInlet(
    const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector <MatrixXd> &ddRindWdW);
#endif
