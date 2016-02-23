#ifndef QUASIONED_H
#define QUASIONED_H

#include<vector>

long double quasiOneD(std::vector <long double> x, 
                 std::vector <long double> dx, 
                 std::vector <long double> S,
                 std::vector <long double> designVar,
                 std::vector <long double> &W);

long double TotalPressureLoss(std::vector <long double> W);

void ioTargetPressure(int io, std::vector <long double> &p);

void inletBC(std::vector <long double> &W, std::vector <long double> &Resi, long double dt0, long double dx0);
void outletBC(std::vector <long double> &W, std::vector <long double> &Resi, long double dt0, long double dx0);

#endif
