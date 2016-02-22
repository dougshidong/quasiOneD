#ifndef FLUX_H
#define FLUX_H

#include<vector>
#include<complex>

void getFlux(std::vector <std::complex<double> > &Flux,
             std::vector <std::complex<double> > W,
             std::vector <std::complex<double> > F);

void Flux_StegerWarming(std::vector <std::complex<double> > &Flux,
                         std::vector <std::complex<double> > W);

void Flux_Scalar(std::vector <std::complex<double> > &Flux,
                 std::vector <std::complex<double> > W,
                 std::vector <std::complex<double> > F);
#endif
