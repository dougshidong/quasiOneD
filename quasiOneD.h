#ifndef QUASIONED_H
#define QUASIONED_H

#include<vector>

double quasiOneD(std::vector <double> x, 
		std::vector <double> dx, 
		std::vector <double> S,
		int fitnessFun,
		std::vector <double> designVar,
		std::vector <double> &W);

double TotalPressureLoss(int nx, std::vector <double> W);

void ioTargetPressure(int io, int nx, std::vector <double> &p);

#endif
