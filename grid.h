#ifndef grid_h
#define grid_h

#include <vector>

void InitializeGrid(int nx, std::vector <double> &x, std::vector <double> &dx,
		std::vector <double> &S);

std::vector <double> calcVolume(std::vector <double> S, std::vector <double> dx);

std::vector <double> evalX(int nx, double a, double b);

std::vector <double> evalDx(int nx, std::vector <double> x);

std::vector <double> evalS(int nx, std::vector <double> geom,
			std::vector <double> x, std::vector <double> dx);



#endif

