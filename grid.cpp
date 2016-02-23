#include <iostream>
#include <math.h>
#include <vector>
#include "globals.h"

// Evaluate X
std::vector <long double> evalX(long double a, long double b)
{
    std::vector <long double> x(nx);
    long double dxConst = (b - a)/nx;
    
    for(int i = 0; i < nx; i++)
        x[i] = dxConst/2 + dxConst * i;

    return x;
}

//  Evaluate dx

std::vector <long double> evalDx(std::vector <long double> x)
{
    std::vector <long double> dx(nx);
    
    dx[0] = x[1] - x[0];
    for(int i = 1; i < nx - 1; i++)
    {
        dx[i] =  (x[i] - x[i - 1])/2  +  (x[i + 1] - x[i])/2 ;
    }
    dx[nx - 1] = x[x.size() - 1] - x[x.size() - 2];

    return dx;
}

std::vector <long double> evalS(std::vector <long double> geom,
                           std::vector <long double> x,
                           std::vector <long double> dx)
{
    std::vector <long double> S(nx + 1);
    long double xh;
    // Define Area
    for(int i = 0; i < nx; i++)
    {
        xh = fabs(x[i] - dx[i] / 2.0);
        S[i] =  1 - geom[0] * pow(sin(PI * pow(xh, geom[1])), geom[2]);
    }
    
    S[nx] =  1 - geom[0] * pow(sin(PI * pow(x[nx - 1] + dx[nx - 1]/2, geom[1])), geom[2]);

    return S;


}

