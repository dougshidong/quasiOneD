#ifndef CONVERT_H
#define CONVERT_H

#include<vector>
#include<Eigen/Core>

void getp(
    std::vector <double> W,
    std::vector <double> &p);

void WtoP(
    std::vector <double> W,
    std::vector <double> &rho,
    std::vector <double> &u,
    std::vector <double> &e);

void WtoP2(
    std::vector <double> W,
    std::vector <double> &rho,
    std::vector <double> &u,
    std::vector <double> &p);

void WtoP(
    std::vector <double> W,
    std::vector <double> &rho,
    std::vector <double> &e,
    std::vector <double> &u,
    std::vector <double> &p,
    std::vector <double> &c,
    std::vector <double> &T);

void dWpdW(
    std::vector <double> &M,
    std::vector <double> W,
    int k);

Eigen::MatrixXd dWpdW(
    std::vector <double> W,
    int i);

void dWdWp(
    std::vector <double> &M,
    std::vector <double> W,
    int k);

void WtoF(
    std::vector <double> W,
    std::vector <double> &F);

void WtoQ(
    std::vector <double> W,
    std::vector <double> &Q,
    std::vector <double> S);
#endif
