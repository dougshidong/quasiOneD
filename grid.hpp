#ifndef GRID_H
#define GRID_H

#include "structures.hpp"
#include <vector>

template<typename dreal>
std::vector<dreal> calcVolume(
    std::vector<dreal> area,
    std::vector<dreal> dx);

template<typename dreal>
std::vector<dreal> uniform_x(dreal a, dreal b, int n_elem);

template<typename dreal>
std::vector<dreal> eval_dx(std::vector<dreal> x);

template<typename dreal>
std::vector<dreal> evalS(
	const struct Design<dreal> &design,
    const std::vector<double> &x,
    const std::vector<double> &dx);

#endif

