#include "quasiOneD.h"
#include "grid.h"
#include "optimizer.h"
#include "adjoint.h"
#include <iostream>
#include <vector>
#include "globals.h"
#include <fenv.h>

int main()
{
    std::vector <long double> x(nx), S(nx + 1);
    std::vector <long double> dx(nx);
    std::vector <long double> geom(3);
    std::vector <long double> W(3 * nx, 0);
    long double fitness;
    feenableexcept(FE_INVALID | FE_OVERFLOW);
    
    geom[0] = h_geom;
    geom[1] = t1_geom;
    geom[2] = t2_geom;
    
    x = evalX(a_geom, b_geom);
    dx = evalDx(x);
    S = evalS(geom, x, dx);

    if(opt == 0)
        fitness = quasiOneD(x, dx, S, geom, W);
    if(opt == 1)
        design(x, dx, S, geom);

    return 0;
}
