#ifndef parametrization_h
#define parametrization_h
#include<vector>
#include<Eigen/Core>
Eigen::MatrixXd evaldAreadDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const struct Design &design);

std::vector <Eigen::MatrixXd> evalddAreadDesdDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const struct Design &design);
#endif

