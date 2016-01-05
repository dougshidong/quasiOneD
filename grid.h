#ifndef grid_h
#define grid_h

#include <vector>

void InitializeGrid(std::vector <double> &x, std::vector <double> &dx,
		std::vector <double> &S);

std::vector <double> calcVolume(std::vector <double> S, std::vector <double> dx);

std::vector <double> evalX(double a, double b);

std::vector <double> evalDx(std::vector <double> x);

std::vector <double> evalS(std::vector <double> geom,
			std::vector <double> x, std::vector <double> dx);



#endif

