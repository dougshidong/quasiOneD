
#ifndef hessianOutlet_h
#define hessianOutlet_h

#include<vector>
#include<Eigen/Core>
#include"structures.h"

using namespace Eigen;

void HessianOutlet(
    const struct Flow_options &flo_opts,
    const std::vector<double> &W,
    std::vector <MatrixXd> &ddRindWdW);
#endif
