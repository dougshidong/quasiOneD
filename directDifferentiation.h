
#ifndef directDifferentiation_h
#define directDifferentiation_h

#include<vector>
#include<Eigen/Core>

MatrixXd evaldWdDes(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> W,
    std::vector <double> designVar);
VectorXd directDifferentiation(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> W,
    std::vector <double> designVar);
#endif
