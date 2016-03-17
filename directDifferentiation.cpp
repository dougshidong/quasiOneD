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
MatrixXd evaldWdDesign(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> W,
    std::vector <double> designVar)
{
    // DR   dR   dR DW           DW   -( dR ) ^ (-1) ( dR )
    // -- = -- + -- -- = 0  -->  -- =  ( -- )        ( -- )
    // DS   dS   dW DS           DS    ( dW )        ( dS )

    //Get Primitive Variables
    std::vector <double> rho(nx), u(nx), e(nx);
    std::vector <double> T(nx), p(nx), c(nx), Mach(nx);
    WtoP(W, rho, u, e, p, c, T);

    // Get Fluxes
    std::vector <double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);

    // Evaluate dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    dRdS = evaldRdS(Flux, S, W);

    // Evaluate dSdDesign
    MatrixXd dSdDesign(nx + 1, designVar.size());
    dSdDesign = evaldSdDesign(x, dx, designVar);

    //Evaluate dRdDesign
    MatrixXd dRdDesign(3 * nx, nDesVar);
    dRdDesign = dRdS * dSdDesign;

    // Evaluate dRdW
    std::vector <double> dt(nx, 1);
    SparseMatrix <double> dRdW;
    dRdW = evaldRdW(W, dx, dt, S, u[0]/c[0]);

    // Solve DWDS
    MatrixXd dWdDesign(3 * nx, nDesVar);
    // Solver type eig_solv
    // 0 = Sparse LU
    // 1 = Dense LU Full Piv
    // 2 = Sparse Iterative BiCGSTAB
    int eig_solv = 0;
    dWdDesign = solveSparseAXB(-dRdW, dRdDesign, eig_solv);

    return dWdDesign;
}

VectorXd directDifferentiation(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> W,
    std::vector <double> designVar)
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
    // Evaluate dSdDesign
    MatrixXd dSdDesign(nx + 1, designVar.size());
    dSdDesign = evaldSdDesign(x, dx, designVar);
    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
    dIcdS.setZero();
    VectorXd dIcdDesign(nDesVar);
    dIcdDesign = dIcdS.transpose() * dSdDesign;

    MatrixXd dWdDesign(3 * nx, nDesVar);
    dWdDesign = evaldWdDesign(x, dx, S, W, designVar);
    // Evaluate DIDS
    VectorXd dIdDesign(nDesVar);
    dIdDesign = (dIcdDesign.transpose() + dIcdW.transpose() * dWdDesign);


    VectorXd grad(designVar.size());
    grad = dIdDesign;
    std::cout<<"Gradient from Direct Differentiation:"<<std::endl;
    std::cout<<std::setprecision(15)<<grad<<std::endl;
    return grad;
}
