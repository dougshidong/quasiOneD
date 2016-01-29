#include <math.h>
#include "globals.h"

const double PI = atan(1.0) * 4.0;

// Number of Cells in Grid
const int nx = 60;

// Optimization Flag
const int opt = 1;

// Number of Design Variables for Optimization
const int nDesVar = 3;
// Input Geometry
const double a_geom = 0, b_geom = 1;
//const double h_geom = 0.025, t1_geom = 0.70, t2_geom = 1.00;
const double h_geom = 0.10, t1_geom = 0.80, t2_geom = 3.00;

// Fitness Function
// 0  -  Total Pressure Loss
// 1  -  Pressure Target
const int fitnessFun = 1;

// Descent Type for Optimization
// 1  -  Steepest Descent
// 4  -  Quasi-Newton (BFGS)
const int descentType = 4;

// Gradient Type for Optimization
// 1  -  FD Forward
// 2  -  FD Backward
// 3  -  FD Centered
// 4  -  Adjoint Method
const int gradientType = 1;

// Create Target Pressure
// 0   -   Do NOT Create Target Pressure
// 1   -   Create Target Pressure
const int createTarget = 0;
