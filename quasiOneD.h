#ifndef QUASIONED_H
#define QUASIONED_H

#include<vector>

double isenP(double pt, double M);
double isenT(double Tt, double M);
double quasiOneD(
    std::vector<double> x,
    std::vector<double> area,
    std::vector<double> &W);
void inletBC(
    std::vector<double> &W,
    std::vector<double> &Resi,
    double dt0, double dx0);
void outletBC(
    std::vector<double> &W,
    std::vector<double> &Resi,
    double dt0, double dx0);

#endif
