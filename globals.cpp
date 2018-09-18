#include <math.h>
#include "globals.h"
#include <string>

const double PI = atan(1.0) * 4.0;

std::string filename;

int n_elem;
double a_geom, b_geom;
double h_geom, t1_geom, t2_geom;

int StepScheme, FluxScheme;
double Scalareps;

double CFL;
double flowConv;
int maxIt;

int printIt, printConv, printW;

double gam, R, Cv;
double a2;
double inlet_mach, inlet_total_T, inlet_total_p, outlet_p;

int opt, desParam, fitnessFun;
int nDesVar;
int descentType, gradientType, hessianType, exactHessian;
double htol;
int nCG;
double gradConv;
int maxDesign;

double h_tar, t1_tar, t2_tar;

int n_control_pts, spline_degree;
