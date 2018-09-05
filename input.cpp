#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "globals.h"
#include "math.h"
#include <string>
#include <algorithm> // remove_if

#define MAX_STRLEN 256
void inputfile()
{
    FILE *inputf;
    inputf = fopen("input.in", "r");

    char buf[MAX_STRLEN];
    // Number of Cells in Grid
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
	printf("asdasd %td %s \n",sizeof(buf), buf);
    filename = buf;
    filename.erase(std::remove_if(filename.begin(), filename.end(), isspace), filename.end());
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%d", &nx);
	printf("asdasd %td %s \n",sizeof(buf), buf);

    // Input Geometry
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%lf %lf", &a_geom, &b_geom);
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%lf %lf %lf", &h_geom, &t1_geom, &t2_geom);

    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
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
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    int tol;
    sscanf(buf, "%lf %d %d", &CFL, &tol, &maxIt);
    flowConv = pow(10.0, tol);

    // Printing Flow Stuff for Debugging
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%d %d %d", &printIt, &printConv, &printW);

    // Flow Inputs
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%lf %lf", &gam, &R);
    Cv = R / (gam - 1.0);

    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    double temp;
    sscanf(buf, "%lf %lf %lf %lf", &Min, &Ttin, &ptin, &temp);
    pexit = temp * ptin;
    a2 = 2.0 * gam * Cv * Ttin * ((gam - 1.0) / (gam + 1.0));

    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    // Design Optimization Parameters
    // opt = 0 or 1 to Turn ON/OFF Optimization
    // Design Variables
    // 0  -  Individual Areas
    // 1  -  Sin Parametrization (Final Project MECH 539)
    // Fitness Function
    // 0  -  Total Pressure Loss
    // 1  -  Pressure Target
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%d %d %d", &opt, &desParam, &fitnessFun);

    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
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
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%d %d %d %d %lf", &descentType, &gradientType, &hessianType, &exactHessian, &htol);

    // Number of CG steps when using Truncated Newton
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%d", &nCG);

    // Design Convergence
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%d %d", &tol, &maxDesign);
    gradConv = pow(10.0, tol);

    // Target Geometry
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%lf %lf %lf", &h_tar, &t1_tar, &t2_tar);

    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Skip Line
    if (fgets(buf, sizeof(buf), inputf) == NULL) {abort();} // Read
    sscanf(buf, "%d %d", &nctl, &spline_degree);

    fclose(inputf);
}

