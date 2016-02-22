#include <math.h>
#include "globals.h"
#include <complex>

const std::complex<double> PI = atan(1.0) * 4.0;

// Number of Cells in Grid
const int nx = 50;

// Optimization Flag
const int opt = 1;

// Number of Design Variables for Optimization
const int nDesVar = 3;
// Input Geometry
const std::complex<double> a_geom = 0, b_geom = 1;
//const std::complex<double> h_geom = 0.10, t1_geom = 0.80, t2_geom = 3.00;
const std::complex<double> h_geom = 0.05, t1_geom = 1.00, t2_geom = 6.00;

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
// 4  -  Complex Step FD Method
const int gradientType = 4;

// Create Target Pressure
// 0   -   Do NOT Create Target Pressure
// 1   -   Create Target Pressure
const int createTarget = 0;


// FLOW RELATED VARIABLES
// Gas Constants
const double gam = 1.4;
const std::complex<double> R = 1716.0 / 1716.0;
const std::complex<double> Cv = R / (gam - 1.0);

// Inlet
const std::complex<double> Min = 0.75;
const std::complex<double> Ttin = 531.2 / 531.2;
const std::complex<double> ptin = 2117.0 / 2117.0;
// Outlet
const std::complex<double> pexit = 0.12 * ptin;
// Constant
const std::complex<double> a2 = 2.0 * gam * Cv * Ttin * ((gam - 1.0)/(gam + 1.0)); // used in isentropic nozzle

// Convergence Settings
const std::complex<double> CFL = 0.2;
const std::complex<double> conv = 1e-14;
const int maxIt = 50000;
const int printIt = 1000;
const int printConv = 0; // 0 to hide real - time convergence
const int printW = 0;

// Stepping Scheme
// 0   -   Euler Explicit
// 1   -   Runge - Kutta 4th order
// 2   -   Jameson's Runge-Kutta 4th order
const int StepScheme = 0;

// Flux Scheme
// 0   -   Steger Warming (SW)
// 1   -   Scalar Dissipation
const int FluxScheme = 1;
std::complex<double> Scalareps = 0.5;
