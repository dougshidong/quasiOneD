#ifndef TIMESTEP_H
#define TIMESTEP_H

#include<vector>

void stepInTime(
    std::vector <double> S,
    std::vector <double> dx,
    std::vector <double> dt,
    std::vector <double> &Resi,
    std::vector <double> &W);

void getDomainResi(
    std::vector <double> W,
    std::vector <double> Flux,
    std::vector <double> S,
    std::vector <double> &Resi);

#endif
