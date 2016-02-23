#ifndef grid_h
#define grid_h

#include <vector>

std::vector <long double> calcVolume(std::vector <long double> S, std::vector <long double> dx);

std::vector <long double> evalX(long double a, long double b);

std::vector <long double> evalDx(std::vector <long double> x);

std::vector <long double> evalS(std::vector <long double> geom,
                           std::vector <long double> x,
                           std::vector <long double> dx);



#endif

