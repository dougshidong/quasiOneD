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
    //
    //  
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);
    ddRdWdS = evalddRdWdS(W, S);
    std::vector <SparseMatrix <double> > ddRdWdSFD(3 * nx);
    ddRdWdSFD = evalddRdWdS_FD(W, S);

    double sum=0;
    for(int i = 4; i < 3 * nx - 3; i++)
    {
        sum += (ddRdWdS[i] - ddRdWdSFD[i]).norm() / (ddRdWdS[i]).norm();
        std::cout<<(ddRdWdS[i] - ddRdWdSFD[i]).norm() / (ddRdWdS[i].norm())<<std::endl;
    }
    std::cout<<"FD vs AN:  "<<sum<<std::endl;


    std::vector < SparseMatrix <double> > ddRdWdW(3 * nx);
    ddRdWdW = evalddRdWdW(W, S);
    
    MatrixXd Hessian(nDesVar, nDesVar);
    return Hessian;
}


