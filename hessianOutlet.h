
#ifndef hessianOutlet_h
#define hessianOutlet_h

#include<vector>
#include<Eigen/Core>

using namespace Eigen;

void HessianOutlet(
    std::vector<double> W,
    std::vector <MatrixXd> &ddRindWdW);
#endif
