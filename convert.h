#ifndef CONVERT_H
#define CONVERT_H

#include<vector>
#include<Eigen/Core>

void getp(
    const std::vector<double> &W,
    std::vector<double> &p);

void WtoP(
    const std::vector<double> &W,
    std::vector<double> &rho,
    std::vector<double> &u,
    std::vector<double> &e);

void WtoP2(
    const std::vector<double> &W,
    std::vector<double> &rho,
    std::vector<double> &u,
    std::vector<double> &p);

void WtoP(
    const std::vector<double> &W,
    std::vector<double> &rho,
    std::vector<double> &e,
    std::vector<double> &u,
    std::vector<double> &p,
    std::vector<double> &c,
    std::vector<double> &T);

void PtoW(
    std::vector<double> &W,
    const std::vector<double> &Wp);

void dWpdW(
    std::vector<double> &M,
    const std::vector<double> &W,
    int k);

Eigen::MatrixXd dWpdW(
    const std::vector<double> &W,
    int i);

void dWdWp(
    std::vector<double> &M,
    const std::vector<double> &W,
    int k);

void WtoF(
    const std::vector<double> &W,
    std::vector<double> &F);

void WtoQ(
    const std::vector<double> &W,
    std::vector<double> &Q,
    const std::vector<double> &S);
#endif
