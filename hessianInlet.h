
#ifndef hessianInlet_h
#define hessianInlet_h

#include<vector>
#include<Eigen/Core>

using namespace Eigen;

void HessianInlet(
    std::vector <double> W,
    std::vector <MatrixXd> &ddRindWdW);
#endif
