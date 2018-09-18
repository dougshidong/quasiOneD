#include <iostream>
#include <vector>
#include <fenv.h>
#include "globals.h"
#include "input.h"
#include "grid.h"
#include "spline.h"
#include "quasiOneD.h"
#include "oneshot.h"
#include "fitness.h"
#include "output.h"
#include "optimizer.h"
#include "convert.h"
#include"petsc.h"
#include"petscsys.h"

static char help[] = "QuasiOneD\n\n";
int main(int argc,char **argv)
{
    inputfile();

    std::vector<double> x(n_elem), area(n_elem + 1);
    std::vector<double> dx(n_elem);
    std::vector<double> W(3 * n_elem, 0);
    feraiseexcept(FE_INVALID | FE_OVERFLOW);

    // Initialize Shape
    std::vector<double> geom(3);

    x = evalX(a_geom, b_geom);
    dx = eval_dx(x);
    if (opt == 0)
    {
        geom[0] = h_geom;
        geom[1] = t1_geom;
        geom[2] = t2_geom;
        area = evalS(geom, x, dx, 1);
        quasiOneD(x, area, W);
    }
    if (opt == 1)
    {
        std::cout<<"Creating Target Pressure"<<std::endl;
        geom[0] = h_tar;
        geom[1] = t1_tar;
        geom[2] = t2_tar;
        area = evalS(geom, x, dx, 1);
        outVec("TargetGeom.dat", "w", x);
        outVec("TargetGeom.dat", "a", area);

        quasiOneD(x, area, W);
        std::vector<double> pt(n_elem);
        getp(W, pt);
        ioTargetPressure(1, pt);

        geom[0] = h_geom;
        geom[1] = t1_geom;
        geom[2] = t2_geom;
        area = evalS(geom, x, dx, 1);
        if (desParam == 0) nDesVar = n_elem - 1;
        if (desParam == 1) nDesVar = 3;
        if (desParam == 2) nDesVar = n_control_pts - 2; // Inlet and Outlet are constant

        std::vector<double> desVar(nDesVar);
        if (desParam == 0)
        {
            for (int Si = 1; Si < n_elem; Si++)
            {
                desVar[Si - 1] = area[Si];
            }
        }
        if (desParam == 1) desVar = geom;
        if (desParam == 2)
        {
            geom = getCtlpts(x, dx, area); 
            for (int iVar = 0; iVar < nDesVar; iVar++)
            {
                desVar[iVar] = geom[iVar+1];
            }
            area = evalS(desVar, x, dx, 2);
        }

        PetscInitialize(&argc, &argv, (char*)0,help);
        design(x, dx, area, desVar);
        PetscFinalize();
    }

    return 0;
}
