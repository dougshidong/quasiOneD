#include <iostream>
#include <math.h>
#include <vector>
#include "globals.h"
#include <complex>

// Evaluate X
std::vector <std::complex<double> > evalX(std::complex<double> a, std::complex<double> b)
{
    std::vector <std::complex<double> > x(nx);
    std::complex<double> dxConst = (b - a)/(double)nx;
    
    for(int i = 0; i < nx; i++)
        x[i] = dxConst/2.0 + dxConst * (double)i;

    return x;
}

//  Evaluate dx

std::vector <std::complex<double> > evalDx(std::vector <std::complex<double> > x)
{
    std::vector <std::complex<double> > dx(nx);
    
    dx[0] = x[1] - x[0];
    for(int i = 1; i < nx - 1; i++)
    {
        dx[i] =  (x[i] - x[i - 1])/2.0  +  (x[i + 1] - x[i])/2.0 ;
    }
    dx[nx - 1] = x[x.size() - 1] - x[x.size() - 2];

    return dx;
}

std::vector <std::complex<double> > evalS(std::vector <std::complex<double> > geom,
                           std::vector <std::complex<double> > x,
                           std::vector <std::complex<double> > dx)
{
    std::vector <std::complex<double> > S(nx + 1);
    std::complex<double> xh;
    // Define Area
    for(int i = 1; i < nx; i++)
    {
        xh = x[i] - dx[i] / 2.0;
        S[i] =  1.0 - geom[0] * pow(sin(PI * pow(xh, geom[1])), geom[2]);
    }
    
    S[0] = 1.0;
    S[nx] =  1.0 - geom[0] * pow(sin(PI * pow(x[nx - 1] + dx[nx - 1]/2.0, geom[1])), geom[2]);

    return S;


}

