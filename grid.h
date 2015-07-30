#ifndef grid_h
#define grid_h

#include <vector>

void InitializeGrid(int nx, std::vector <double> &x, std::vector <double> &dx,
		std::vector <double> &S);

std::vector <double> calcVolume(std::vector <double> S, std::vector <double> dx);


#endif

