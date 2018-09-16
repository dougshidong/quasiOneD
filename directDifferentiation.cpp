#include<math.h>
#include<vector>
#include<iostream>
#include<iomanip>
#include<Eigen/Eigen>
#include "globals.h"
#include "adjoint.h"
#include "convert.h"
#include "flux.h"
#include "residuald1.h"
#include "parametrization.h"


using namespace Eigen;
MatrixXd evaldWdDes(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> S,
    std::vector<double> W,
    std::vector<double> designVar)
{
    // DR   dR   dR DW           DW   -( dR ) ^ (-1) ( dR )
    // -- = -- + -- -- = 0  -->  -- =  ( -- )        ( -- )
    // DS   dS   dW DS           DS    ( dW )        ( dS )

    //Get Primitive Variables
    std::vector<double> rho(nx), u(nx), e(nx);
    std::vector<double> T(nx), p(nx), c(nx), Mach(nx);
    WtoP(W, rho, u, e, p, c, T);

    // Get Fluxes
    std::vector<double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);

    // Evaluate dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    dRdS = evaldRdS(Flux, S, W);

    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, designVar.size());
    dSdDes = evaldSdDes(x, dx, designVar);

    //Evaluate dRdDes
    MatrixXd dRdDes(3 * nx, nDesVar);
    dRdDes = dRdS * dSdDes;

    // Evaluate dRdW
    std::vector<double> dt(nx, 1);
    SparseMatrix<double> dRdW;
    dRdW = evaldRdW(W, dx, dt, S);

    // Solve DWDS
    MatrixXd dWdDes(3 * nx, nDesVar);
    // Solver type eig_solv
    // 0 = Sparse LU
    // 1 = Dense LU Full Piv
    // 2 = Sparse Iterative BiCGSTAB
    int eig_solv = 0;
    dWdDes = solveSparseAXB(-dRdW, dRdDes, eig_solv);

    return dWdDes;
}

VectorXd directDifferentiation(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> S,
    std::vector<double> W,
    std::vector<double> designVar)
{
    // Direct Differentiation
    // I = Ic(W, S)
    // R = R(W, S) = 0 @ SS
    // W = W(W0, S)
    // DI   dIc   dIc DW
    // -- = --- + --- --
    // DS   dS    dW  DS
    //
    // Evaluate dIcdW
    VectorXd dIcdW(3 * nx);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, designVar.size());
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
    dIcdS.setZero();
    VectorXd dIcdDes(nDesVar);
    dIcdDes = dIcdS.transpose() * dSdDes;

    // Evaluate dWdDes
    MatrixXd dWdDes(3 * nx, nDesVar);
    dWdDes = evaldWdDes(x, dx, S, W, designVar);
    // Evaluate dIcdDes
    VectorXd dIdDes(nDesVar);
    dIdDes = (dIcdDes.transpose() + dIcdW.transpose() * dWdDes);


    VectorXd grad(designVar.size());
    grad = dIdDes;
    std::cout<<"Gradient from Direct Differentiation:"<<std::endl;
    std::cout<<std::setprecision(15)<<grad<<std::endl;
    return grad;
}
