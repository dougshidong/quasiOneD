#ifndef FLUX_H
#define FLUX_H

#include<vector>

void getFlux(std::vector <long double> &Flux,
             std::vector <long double> W,
             std::vector <long double> F);

void Flux_StegerWarming(std::vector <long double> &Flux,
                         std::vector <long double> W);

void Flux_Scalar(std::vector <long double> &Flux,
                 std::vector <long double> W,
                 std::vector <long double> F);
#endif
