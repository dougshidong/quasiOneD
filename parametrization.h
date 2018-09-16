#ifndef parametrization_h
#define parametrization_h
#include<Eigen/Core>
Eigen::MatrixXd evaldSdDes(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> designVar);
std::vector <Eigen::MatrixXd> evalddSdDesdDes(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> designVar);
#endif

