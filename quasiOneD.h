#ifndef QUASIONED_H
#define QUASIONED_H

#include<vector>

double quasiOneD(int nx, 
		std::vector <double> x, 
		std::vector <double> dx, 
		std::vector <double> S,
		int fitnessFun);

double TotalPressureLoss(int nx, std::vector <double> W);

void ioTargetPressure(int io, int nx, std::vector <double> &p);

#endif
