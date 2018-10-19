#ifndef PARAMETRIZATION_H
#define PARAMETRIZATION_H
#include"structures.hpp"
#include<vector>
#include<Eigen/Core>

Eigen::MatrixXd evaldAreadDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const struct Design<double> &design);

Eigen::MatrixXd evaldAreadDes_FD(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const struct Design<double> &design);

std::vector <Eigen::MatrixXd> evalddAreadDesdDes(
    const std::vector<double> &x,
    const std::vector<double> &dx,
    const struct Design<double> &design);
#endif

