#ifndef grid_h
#define grid_h

#include <vector>
#include<complex>

std::vector <std::complex<double> > calcVolume(std::vector <std::complex<double> > S, std::vector <std::complex<double> > dx);

std::vector <std::complex<double> > evalX(std::complex<double> a, std::complex<double> b);

std::vector <std::complex<double> > evalDx(std::vector <std::complex<double> > x);

std::vector <std::complex<double> > evalS(std::vector <std::complex<double> > geom,
                           std::vector <std::complex<double> > x,
                           std::vector <std::complex<double> > dx);



#endif

