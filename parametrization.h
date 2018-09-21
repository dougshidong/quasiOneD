#ifndef parametrization_h
#define parametrization_h
#include<vector>
#include<Eigen/Core>
Eigen::MatrixXd evaldSdDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const struct Design &design);

std::vector <Eigen::MatrixXd> evalddSdDesdDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const std::vector<double> &designVar);
#endif

