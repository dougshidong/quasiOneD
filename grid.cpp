#include <iostream>
#include <math.h>
#include <vector>
#include "globals.h"
#include "spline.h"

std::vector<double> sinParam(
    std::vector<double> geom,
    std::vector<double> x,
    std::vector<double> dx);


// Evaluate X
std::vector<double> evalX(double a, double b) {
    std::vector<double> x(nx);
    double dxConst = (b - a)/nx;

    for (int i = 0; i < nx; i++)
        x[i] = dxConst/2 + dxConst * i;

    return x;
}

//  Evaluate dx

std::vector<double> evalDx(std::vector<double> x) {
    std::vector<double> dx(nx);

    dx[0] = x[1] - x[0];
    for (int i = 1; i < nx - 1; i++) {
        dx[i] =  (x[i] - x[i - 1])/2  +  (x[i + 1] - x[i])/2 ;
    }
    dx[nx - 1] = x[x.size() - 1] - x[x.size() - 2];

    return dx;
}

std::vector<double> evalS(
    std::vector<double> geom,
    std::vector<double> x,
    std::vector<double> dx,
    int param)
{
    std::vector<double> area(nx + 1);
    if (param == 0) {
        area[0] = 1.0;
        area[nx] = 1.0;
        for (int Si = 1; Si < nx; Si++)
        {
            area[Si] = geom[Si - 1];
        }
    }
    else if (param == 1)
    {
        area = sinParam(geom, x, dx);\
    }
    else if (param == 2)
    {
        std::vector<double> control_pts(n_control_pts);
        for (int ictl = 1; ictl < n_control_pts - 1; ictl++)
        {
            control_pts[ictl] = geom[ictl - 1];
        }
        control_pts[0] = 1;
        control_pts[n_control_pts - 1] = 1;
        area = evalSpline(control_pts, x, dx);
    }
    return area;
}

std::vector<double> sinParam(
    std::vector<double> geom,
    std::vector<double> x,
    std::vector<double> dx)
{
    std::vector<double> area(nx + 1);
    double xh;
    // Define Area
    for (int i = 1; i < nx; i++)
    {
        xh = x[i] - dx[i] / 2.0;
        area[i] = 1 - geom[0] * pow(sin(PI * pow(xh, geom[1])), geom[2]);
    }

    area[0] = 1;
    area[nx] = 1;

    return area;
}

