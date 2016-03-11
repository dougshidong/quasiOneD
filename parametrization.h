#ifndef parametrization_h
#define parametrization_h
#include<Eigen/Core>
Eigen::MatrixXd evaldSdDesign(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> designVar);
#endif

