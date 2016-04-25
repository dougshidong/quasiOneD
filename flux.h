#ifndef FLUX_H
#define FLUX_H

#include<vector>

void getFlux(
    std::vector <double> &Flux,
    const std::vector <double> &W);

void initializeFlux(int nx);

void Flux_Scalar(
    std::vector <double> &Flux,
    const std::vector <double> &W);

void Flux_SW(
    std::vector <double> &Flux,
    const std::vector <double> &W);

void Flux_MSW(
    std::vector <double> &Flux,
    const std::vector <double> &W);

void Flux_CMSW(
    std::vector <double> &Flux,
    const std::vector <double> &W);

void Flux_Roe(
    std::vector <double> &Flux,
    const std::vector <double> &W);
#endif
