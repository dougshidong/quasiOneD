#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>

template<class T> void UNUSED( const T& ) { }

extern const double PI;

extern std::string filename;

extern int nx;
extern double a_geom, b_geom;
extern double h_geom, t1_geom, t2_geom;

extern int StepScheme, FluxScheme;
extern double Scalareps;

extern double CFL;
extern double flowConv;
extern int maxIt;

extern int printIt, printConv, printW;

extern double gam, R, Cv;
extern double Min, Ttin, ptin, pexit;
extern double a2;

extern int opt, desParam, fitnessFun;
extern int nDesVar;
extern int descentType, gradientType, hessianType, exactHessian;
extern double htol;
extern int nCG;
extern double gradConv;
extern int maxDesign;

extern double h_tar, t1_tar, t2_tar;

extern int nctl, spline_degree;

#endif
