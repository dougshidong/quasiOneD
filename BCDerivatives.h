#ifndef BCDerivative_h
#define BCDerivative_h

#include<vector>
#include<Eigen/Sparse>

using namespace Eigen;

void HessianBC_FD(
    std::vector <double> W,
    std::vector <MatrixXd> &ddRindWdW,
    std::vector <MatrixXd> &ddRoutdWdW);

void BCJac(
    std::vector <double> W,
    std::vector <double> dt,
    std::vector <double> dx,
    std::vector <double> &dBidWi,
    std::vector <double> &dBidWd,
    std::vector <double> &dBodWd,
    std::vector <double> &dBodWo);

#endif
