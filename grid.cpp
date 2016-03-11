#include <iostream>
#include <math.h>
#include <vector>
#include "globals.h"

std::vector <double> sinParam(
    std::vector <double> geom,
    std::vector <double> x,
    std::vector <double> dx);

// Evaluate X
std::vector <double> evalX(double a, double b)
{
    std::vector <double> x(nx);
    double dxConst = (b - a)/nx;

    for(int i = 0; i < nx; i++)
        x[i] = dxConst/2 + dxConst * i;

    return x;
}

//  Evaluate dx

std::vector <double> evalDx(std::vector <double> x)
{
    std::vector <double> dx(nx);

    dx[0] = x[1] - x[0];
    for(int i = 1; i < nx - 1; i++)
    {
        dx[i] =  (x[i] - x[i - 1])/2  +  (x[i + 1] - x[i])/2 ;
    }
    dx[nx - 1] = x[x.size() - 1] - x[x.size() - 2];

    return dx;
}

std::vector <double> evalS(
    std::vector <double> geom,
    std::vector <double> x,
    std::vector <double> dx,
    int param)
{
    std::vector <double> S(nx + 1);
    if(param == 0) S = geom;
    else if(param == 1) S = sinParam(geom, x, dx);
    return S;
}

std::vector <double> sinParam(
    std::vector <double> geom,
    std::vector <double> x,
    std::vector <double> dx)
{
    std::vector <double> S(nx + 1);
    double xh;
    // Define Area
    for(int i = 1; i < nx; i++)
    {
        xh = x[i] - dx[i] / 2.0;
        S[i] = 1 - geom[0] * pow(sin(PI * pow(xh, geom[1])), geom[2]);
    }

    S[0] = 1;
    S[nx] = 1;

    return S;
}
