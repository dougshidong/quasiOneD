#ifndef QUASIONED_H
#define QUASIONED_H

#include<vector>
#include <complex>

std::complex<double> quasiOneD(std::vector <std::complex<double> > x, 
                 std::vector <std::complex<double> > dx, 
                 std::vector <std::complex<double> > S,
                 std::vector <std::complex<double> > designVar,
                 std::vector <std::complex<double> > &W);

std::complex<double> TotalPressureLoss(std::vector <std::complex<double> > W);

void ioTargetPressure(int io, std::vector <std::complex<double> > &p);

void inletBC(std::vector <std::complex<double> > &W, std::vector <std::complex<double> > &Resi, std::complex<double> dt0, std::complex<double> dx0);
void outletBC(std::vector <std::complex<double> > &W, std::vector <std::complex<double> > &Resi, std::complex<double> dt0, std::complex<double> dx0);

#endif
