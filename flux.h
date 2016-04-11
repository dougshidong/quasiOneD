#ifndef FLUX_H
#define FLUX_H

#include<vector>

void getFlux(
    std::vector <double> &Flux,
    std::vector <double> W);

void Flux_StegerWarming(
    std::vector <double> &Flux,
    std::vector <double> W);

void Flux_Scalar(
    std::vector <double> &Flux,
    std::vector <double> W);

void Flux_SW(
    std::vector <double> &Flux,
    std::vector <double> W);
#endif
