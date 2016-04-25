#ifndef TIMESTEP_H
#define TIMESTEP_H

#include<vector>

void initializeTimeStep(int nx);

void stepInTime(
    const std::vector <double> &S,
    const std::vector <double> &dx,
    const std::vector <double> &dt,
    std::vector <double> &Resi,
    std::vector <double> &W);

void getDomainResi(
    const std::vector <double> &W,
    const std::vector <double> &Flux,
    const std::vector <double> &S,
    std::vector <double> &Resi);

#endif
