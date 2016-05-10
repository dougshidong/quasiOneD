#include <iostream>
#include <stdio.h>
#include "globals.h"
#include "math.h"
#include <string>
#include <algorithm> // remove_if

void inputfile()
{
    FILE *inputf;
    inputf = fopen("input.in", "r");

    char buf[100];
    // Number of Cells in Grid
    fgets(buf, sizeof buf, inputf); // Skip Line
    filename = buf;
    filename.erase(std::remove_if(filename.begin(), filename.end(), isspace), filename.end());
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%d", &nx);

    // Input Geometry
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%lf %lf", &a_geom, &b_geom);
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%lf %lf %lf", &h_geom, &t1_geom, &t2_geom);

    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    // Flow Solver Paramameter
    // Stepping Scheme
    // 0   -   Euler Explicit
    // 1   -   Runge - Kutta 4th order
    // 2   -   Jameson's Runge-Kutta 4th order
    // Flux Scheme
    // 0   -   Steger Warming (SW)
    // 1   -   Scalar Dissipation (SD)
    // Scalareps only necessary for SD
    sscanf(buf, "%d %d %lf", &StepScheme, &FluxScheme, &Scalareps);

    // Flow Convergence
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    int tol;
    sscanf(buf, "%lf %d %d", &CFL, &tol, &maxIt);
    flowConv = pow(10.0, tol);

    // Printing Flow Stuff for Debugging
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%d %d %d", &printIt, &printConv, &printW);

    // Flow Inputs
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%lf %lf", &gam, &R);
    Cv = R / (gam - 1.0);

    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    double temp;
    sscanf(buf, "%lf %lf %lf %lf", &Min, &Ttin, &ptin, &temp);
    pexit = temp * ptin;
    a2 = 2.0 * gam * Cv * Ttin * ((gam - 1.0) / (gam + 1.0));

    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    // Design Optimization Parameters
    // opt = 0 or 1 to Turn ON/OFF Optimization
    // Design Variables
    // 0  -  Individual Areas
    // 1  -  Sin Parametrization (Final Project MECH 539)
    // Fitness Function
    // 0  -  Total Pressure Loss
    // 1  -  Pressure Target
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%d %d %d", &opt, &desParam, &fitnessFun);

    fgets(buf, sizeof buf, inputf); // Skip Line
    // Descent Type for Optimization
    // 1  -  Steepest Descent
    // 2  -  Quasi-Newton (BFGS)
    // 3  -  Newton
    // 4  -  Truncated-Newton
    // Gradient Type
    //-3  -  FD Centered
    //-2  -  FD Centered
    //-1  -  FD Centered
    // 1  -  Adjoint Variable
    // 2  -  Direct Differentiation
    // Hessian Type
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%d %d %d %d", &descentType, &gradientType, &hessianType, &exactHessian);

    // Number of CG steps when using Truncated Newton
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%d", &nCG);

    // Design Convergence
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%d %d", &tol, &maxDesign);
    gradConv = pow(10.0, tol);

    // Target Geometry
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%lf %lf %lf", &h_tar, &t1_tar, &t2_tar);

    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Skip Line
    fgets(buf, sizeof buf, inputf); // Read
    sscanf(buf, "%d %d", &nctl, &spline_degree);

    fclose(inputf);
}

