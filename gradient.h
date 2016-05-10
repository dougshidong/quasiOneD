#ifndef gradient_h
#define gradient_h
#include <Eigen/Core>
#include <vector>
VectorXd getGradient(int gType,
    double currentI,
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> W,
    std::vector <double> designVar,
    VectorXd &psi);
#endif
