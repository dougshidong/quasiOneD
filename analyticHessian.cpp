#include<iostream>
#include<math.h>
#include<vector>
#include<Eigen/Eigen>
#include "globals.h"
#include "residuald2.h"

using namespace Eigen;

MatrixXd directDirectHessian(
    std::vector <double> W,
    std::vector <double> S);

MatrixXd getAnalyticHessian(
    std::vector <double> W,
    std::vector <double> S,
    int method)
{
    MatrixXd Hessian(nDesVar, nDesVar);

    Hessian = directDirectHessian(W, S);

    return Hessian;
}

MatrixXd directDirectHessian(
    std::vector <double> W,
    std::vector <double> S)
{
    // I = Ic(W, S)
    // R = R(W, S) = 0
    // W = W(S)
    //
    // 1st Derivative
    // DI    dIc   dIc dW
    // --  = --- + --- ---
    // DSi   dSi   dW  dSi
    //
    // 2nd Derivative
    //  DDI      ddIc    ddIc  dW    ( dIc    ddIc dW ) dW    dIc ( ddW      ddW  dW )
    // ------ = ------ + ----- --- + (----- + ---- ---) --- + --- (------ + ----- ---)
    // DSiDSj   dSidSj   dSidW dSj   (dWdSj   dWdW dSj) dSi   dW  (dSidSj   dSidW dSj)
    //

    //Get Primitive Variables
    std::vector <double> rho(nx), u(nx), e(nx);
    std::vector <double> T(nx), p(nx), c(nx), Mach(nx);
    WtoP(W, rho, u, e, p, c, T);

    // Evaluate dRdW
    std::vector <double> dt(nx, 1);
    SparseMatrix <double> dRdW;
    dRdW = evaldRdW(W, dx, dt, S, u[0]/c[0]);

    // Evaluate ddRdWdS
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);
    ddRdWdS = evalddRdWdS(W, S);

    // Evaluate ddRdWdW
    std::vector < SparseMatrix <double> > ddRdWdW(3 * nx);
    ddRdWdW = evalddRdWdW(W, S);

//  std::vector < SparseMatrix <double> > ddRdWdW_FD(3 * nx);
//  ddRdWdW_FD = evalddRdWdW_FD(W, S);
//  Confirmed domain ddRdWdW through Finite Difference
//  for(int i = 4; i < 3 * nx - 3; i++)
//  {
//      std::cout<<(ddRdWdW[i] - ddRdWdW_FD[i]).norm() / (ddRdWdW[i].norm())<<std::endl;
//  }
//  int outi = 29, outk = 1;
//  MatrixXd matout1 = ddRdWdW[outi*3 + outk].block((outi-1)*3,(outi-1)*3,9,9);
//  std::cout<<matout1<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<std::endl;
//  MatrixXd matout2 = ddRdWdW_FD[outi*3 + outk].block((outi-1)*3,(outi-1)*3,9,9);
//  std::cout<<matout2<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<std::endl;
//  std::cout<<std::endl;
//  MatrixXd matout3 = (matout1 - matout2);
//  std::cout<<matout3<<std::endl;


    MatrixXd Hessian(nDesVar, nDesVar);
    return Hessian;
}


