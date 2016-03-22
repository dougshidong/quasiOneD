#ifndef QUASIONED_H
#define QUASIONED_H

#include<vector>

double quasiOneD(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> &W);

double TotalPressureLoss(std::vector <double> W);

void ioTargetPressure(int io, std::vector <double> &p);

void inletBC(
    std::vector <double> &W,
    std::vector <double> &Resi,
    double dt0, double dx0);
void outletBC(
    std::vector <double> &W,
    std::vector <double> &Resi,
    double dt0, double dx0);

#endif
