#ifndef grid_h
#define grid_h

#include "structures.h"
#include <vector>

std::vector<double> calcVolume(
    std::vector<double> area,
    std::vector<double> dx);

std::vector<double> uniform_x(double a, double b, int n_elem);

std::vector<double> eval_dx(std::vector<double> x);

std::vector<double> evalS(
	const struct Design &design,
    const std::vector<double> &x,
    const std::vector<double> &dx);



#endif

