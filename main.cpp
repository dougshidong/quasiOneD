#include "quasiOneD.h"
#include "grid.h"
#include "optimizer.h"
#include "adjoint.h"
#include <iostream>
#include <stdio.h>
#include <vector>
#include "globals.h"
#include <fenv.h>

int main()
{
    std::vector <double> x(nx), S(nx + 1);
    std::vector <double> dx(nx);
    std::vector <double> W(3 * nx, 0);
    double fitness;
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    // Initialize Shape
    std::vector <double> geom(3);
    geom[0] = h_geom;
    geom[1] = t1_geom;
    geom[2] = t2_geom;
    x = evalX(a_geom, b_geom);
    dx = evalDx(x);
    S = evalS(geom, x, dx, 1);

    if(opt == 0) fitness = quasiOneD(x, dx, S, W);
    if(opt == 1)
    {
        if(desParam == 0) nDesVar = nx + 1;
        if(desParam == 1) nDesVar = 3;
        std::vector <double> desVar(nDesVar);
        if(desParam == 0) desVar = S;
        if(desParam == 1) desVar = geom;
        design(x, dx, S, desVar);
    }

    int nI = 25;
    double hs = 0.04, he = 0.11, dh = (he - hs) / (nI - 1);
    double t1s = 0.7, t1e = 1.1, dt1 = (t1e - t1s) / (nI - 1);
    double t2s = 1.0, t2e = 9.0, dt2 = (t2e - t2s) / (nI - 1);

//  FILE  * Results;
//  Results = fopen("CostPlot.dat", "w");
//  fprintf(Results, "%d\n", nI);
//  for(int i = 0; i < nI; i++)
//      fprintf(Results, "%.15f\n", hs + i * dh);
//  for(int i = 0; i < nI; i++)
//      fprintf(Results, "%.15f\n", t1s + i * dt1);
//  for(int i = 0; i < nI; i++)
//      fprintf(Results, "%.15f\n", t2s + i * dt2);

//  for(int i = 0; i < nI; i++)
//  {
//      std::cout<<"i: "<<i<<std::endl;
//      for(int j = 0; j < nI; j++)
//      {
//          for(int k = 0; k < nI; k++)
//          {
//              std::cout<<"j: "<<j<<std::endl;
//              geom[0] = hs + i * dh;
//              geom[1] = t1s + j * dt1;
//              geom[2] = t2_geom + k * dt2;
//              S = evalS(geom, x, dx);
//              
//              fitness = quasiOneD(x, dx, S, W);
//              fprintf(Results, "%.15f\n", fitness);
//          }
//      }
//  }

//  fclose(Results);
    return 0;
}
