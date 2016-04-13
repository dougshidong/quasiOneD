#include <iostream>
#include <stdio.h>
#include <vector>
#include <fenv.h>
#include "globals.h"
#include "input.h"
#include "grid.h"
#include "quasiOneD.h"
#include "optimizer.h"
#include "convert.h"
#include"petsc.h"
#include"petscsys.h"

static char help[] = "QuasiOneD\n\n";
int main(int argc,char **argv)
{
    inputfile();

    std::vector <double> x(nx), S(nx + 1);
    std::vector <double> dx(nx);
    std::vector <double> W(3 * nx, 0);
    double fitness;
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    // Initialize Shape
    std::vector <double> geom(3);

    x = evalX(a_geom, b_geom);
    dx = evalDx(x);
    if(opt == 0)
    {
        geom[0] = h_geom;
        geom[1] = t1_geom;
        geom[2] = t2_geom;
        S = evalS(geom, x, dx, 1);
        fitness = quasiOneD(x, dx, S, W);
    }
    if(opt == 1)
    {
        std::cout<<"Creating Target Pressure"<<std::endl;
        geom[0] = h_tar;
        geom[1] = t1_tar;
        geom[2] = t2_tar;
        S = evalS(geom, x, dx, 1);

        fitness = quasiOneD(x, dx, S, W);
        std::vector <double> pt(nx);
        getp(W, pt);
        ioTargetPressure(1, pt);

        geom[0] = h_geom;
        geom[1] = t1_geom;
        geom[2] = t2_geom;
        S = evalS(geom, x, dx, 1);
        if(desParam == 0) nDesVar = nx - 1;
        if(desParam == 1) nDesVar = 3;

        std::vector <double> desVar(nDesVar);
        if(desParam == 0)
        {
            for(int Si = 1; Si < nx; Si++)
            {
                desVar[Si - 1] = S[Si];
            }
        }

        if(desParam == 1) desVar = geom;

        PetscInitialize(&argc, &argv, (char*)0,help);
        design(x, dx, S, desVar);
        PetscFinalize();
    }

    return 0;
}
