#include<iostream>
#include<iomanip>
#include<math.h>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include"globals.h"
#include"convert.h"
#include"parametrization.h"
#include"residuald1.h"
#include"residuald2.h"
#include"objectiveDerivatives.h"
#include"BCDerivatives.h"
#include"directDifferentiation.h"
#include"flux.h"
#include"quasiOneD.h"
#include"grid.h"
#include"petscGMRES.h"

using namespace Eigen;

std::vector <MatrixXd> evalddWdDesdDes(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar);

std::vector <MatrixXd> evalddWdDesdDes_FD(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> area,
    std::vector<double> designVar);

MatrixXd directAdjointHessian(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar);
MatrixXd adjointAdjointHessian(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar);
MatrixXd adjointDirectHessian(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar);
MatrixXd directDirectHessian(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar);

MatrixXd getAnalyticHessian(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar,
    int method)
{
    MatrixXd Hessian(nDesVar, nDesVar);

    switch(method) {
        case 0: Hessian = directDirectHessian(x, dx, W, area, designVar);
        case 1: Hessian = adjointDirectHessian(x, dx, W, area, designVar);
        case 2: Hessian = adjointAdjointHessian(x, dx, W, area, designVar);
        case 3: Hessian = directAdjointHessian(x, dx, W, area, designVar);
    }

    return Hessian;
}

MatrixXd directAdjointHessian(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(n_elem + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(n_elem + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dIcdW
    VectorXd dIcdW(3 * n_elem);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate ddIcdWdW
    SparseMatrix<double> ddIcdWdW;
    ddIcdWdW = evaldIcdWdW(W, dx);
    // Evaluate ddIcdWdDes
    MatrixXd ddIcdWdDes(3 * n_elem, nDesVar);
    ddIcdWdDes = evalddIcdWdS() * dSdDes;
    // Evaluate dIcdS
    VectorXd dIcdS(n_elem + 1);
    dIcdS = evaldIcdS();
    // Evaluate dIcdDes
    VectorXd dIcdDes(nDesVar);
    dIcdDes = dIcdS.transpose() * dSdDes;
    // Evaluate ddIcdDesdDes
    MatrixXd ddIcdDesdDes(nDesVar, nDesVar);
    ddIcdDesdDes = dSdDes.transpose() * evalddIcdSdS() * dSdDes;
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    //// Get Fluxes
    std::vector<double> Flux(3 * (n_elem + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * n_elem, n_elem + 1);
    dRdS = evaldRdS(Flux, area, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * n_elem, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector<double> dt(n_elem, 1);
    SparseMatrix<double> dRdW;
    dRdW = evaldRdW(W, dx, dt, area);
    // Evaluate ddRdWdS
    std::vector <SparseMatrix<double> > ddRdWdS(3 * n_elem);// by (3 * n_elem, n_elem + 1)
    ddRdWdS = evalddRdWdS(W, area);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * n_elem);// by (3 * n_elem, nDesVar)
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
       ddRdWdDes[Ri] = ddRdWdS[Ri] * dSdDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(W, area);

    // *************************************
    // Solve for Adjoint (1 Flow Eval)
    // *************************************
    VectorXd psi(3 * n_elem);
      SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
      slusolver1.compute(-dRdW.transpose());
      if (slusolver1.info() != 0)
          std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
    psi = slusolver1.solve(dIcdW);
    //psi = solveGMRES(-dRdW.transpose(),dIcdW);

    // *************************************
    // Evaluate dWdDes (nDesVar Flow Eval)
    // *************************************
    MatrixXd dWdDes(3 * n_elem, nDesVar);
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > factdrdw;
    factdrdw.compute(-dRdW);
    if (factdrdw.info() != 0)
        std::cout<<"Factorization failed. Error: "<<factdrdw.info()<<std::endl;
    if (exactHessian == 1)  dWdDes = factdrdw.solve(dRdDes);
    else if (exactHessian < 0)
    {
        // Iterative Solution of dWdDes
        dWdDes = solveGMRES(-dRdW,dRdDes);

        std::cout<<"dWdDes ||r||/||b|| residual:"<<std::endl;
        std::cout<<(-dRdW*dWdDes - dRdDes).norm()/dRdDes.norm()<<std::endl;

        // Direct Solution of dWdDes
        MatrixXd realdWdDes(3*n_elem,nDesVar);
        realdWdDes = factdrdw.solve(dRdDes);

        std::cout<<"Relative error of approximate dWdDes vs exact dWdDes:"<<std::endl;
        std::cout<<(realdWdDes - dWdDes).norm()/realdWdDes.norm()<<std::endl;
        for (int icol = 0; icol < nDesVar; icol++) {
            std::cout << icol << "\t" << (realdWdDes.col(icol) - dWdDes.col(icol)).norm()/realdWdDes.col(icol).norm()<<std::endl;
        }
    }


    // *************************************
    // Evaluate total derivative DDIcDDesDDes
    // *************************************
    MatrixXd DDIcDDesDDes(nDesVar, nDesVar);
    DDIcDDesDDes.setZero();
    DDIcDDesDDes = ddIcdDesdDes;
    DDIcDDesDDes += dWdDes.transpose() * ddIcdWdDes;
    DDIcDDesDDes += (dWdDes.transpose() * ddIcdWdDes).transpose();
    DDIcDDesDDes += dWdDes.transpose() * ddIcdWdW * dWdDes;
    for (int Si = 0; Si < n_elem + 1; Si++) {
        DDIcDDesDDes += dIcdS(Si) * ddSdDesdDes[Si];
        DDIcDDesDDes += psi.dot(dRdS.col(Si)) * ddSdDesdDes[Si];
    }
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
//        DDIcDDesDDes += psi(Ri) * ddRdDesdDes; //ddRdDesdDes is 0
        DDIcDDesDDes += psi(Ri) * (dWdDes.transpose() * ddRdWdDes[Ri]);
        DDIcDDesDDes += (psi(Ri) * (dWdDes.transpose() * ddRdWdDes[Ri])).transpose();
        DDIcDDesDDes += psi(Ri) * (dWdDes.transpose() * ddRdWdW[Ri] * dWdDes);
    }

    return DDIcDDesDDes;
}
MatrixXd adjointAdjointHessian(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(n_elem + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(n_elem + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dIcdW
    VectorXd dIcdW(3 * n_elem);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate ddIcdWdW
    SparseMatrix<double> ddIcdWdW;
    ddIcdWdW = evaldIcdWdW(W, dx);
    // Evaluate ddIcdWdDes
    MatrixXd ddIcdWdDes(3 * n_elem, nDesVar);
    ddIcdWdDes = evalddIcdWdS() * dSdDes;
    // Evaluate dIcdS
    VectorXd dIcdS(n_elem + 1);
    dIcdS = evaldIcdS();
    // Evaluate dIcdDes
    VectorXd dIcdDes(nDesVar);
    dIcdDes = dIcdS.transpose() * dSdDes;
    // Evaluate ddIcdDesdDes
    MatrixXd ddIcdDesdDes(nDesVar, nDesVar);
    ddIcdDesdDes = dSdDes.transpose() * evalddIcdSdS() * dSdDes;
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    //// Get Fluxes
    std::vector<double> Flux(3 * (n_elem + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * n_elem, n_elem + 1);
    dRdS = evaldRdS(Flux, area, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * n_elem, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector<double> dt(n_elem, 1);
    SparseMatrix<double> dRdW;
    dRdW = evaldRdW(W, dx, dt, area);
    // Evaluate ddRdWdS
    std::vector <SparseMatrix<double> > ddRdWdS(3 * n_elem);// by (3 * n_elem, n_elem + 1)
    ddRdWdS = evalddRdWdS(W, area);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * n_elem);// by (3 * n_elem, nDesVar)
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
       ddRdWdDes[Ri] = ddRdWdS[Ri] * dSdDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(W, area);

    // *************************************
    // Sparse LU of Jacobian Transpose dRdW
    // *************************************
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
    slusolver1.compute(-dRdW.transpose());
    if (slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;

    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver2;
    slusolver2.compute(-dRdW);
    if (slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver2.info()<<std::endl;
    // *************************************
    // Solve for Adjoint 1 psi(1 Flow Eval)
    // *************************************
    VectorXd psi(3 * n_elem);
    psi = slusolver1.solve(dIcdW);

    // *************************************
    // Solve for Adjoint 2 lambda (nDesVar Flow Eval)
    // *************************************
    MatrixXd lambda(3 * n_elem, nDesVar);
    lambda = slusolver2.solve(dRdDes);

    // *************************************
    // Solve for Adjoint 3 eta (nDesVar Flow Eval)
    // *************************************
    MatrixXd eta(3 * n_elem, nDesVar);
    MatrixXd RHS(3 * n_elem, nDesVar);

    RHS = ddIcdWdDes + ddIcdWdW * lambda;
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        RHS += psi(Ri) * ddRdWdDes[Ri];
        RHS += psi(Ri) * ddRdWdW[Ri] * lambda;
    }
    eta = slusolver1.solve(RHS);

    // *************************************
    // Evaluate total derivative DDIcDDesDDes
    // *************************************
    MatrixXd DDIcDDesDDes(nDesVar, nDesVar);
    DDIcDDesDDes.setZero();
    DDIcDDesDDes = ddIcdDesdDes
                   + lambda.transpose() * ddIcdWdDes
                   + eta.transpose() * dRdDes;
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
//      DDIcDDesDDes += psi(Ri) * ddRdDesdDes[Ri];
        DDIcDDesDDes += psi(Ri) * (lambda.transpose() * ddRdWdDes[Ri]);
    }
    for (int Si = 0; Si < n_elem + 1; Si++) {
        DDIcDDesDDes += dIcdS(Si) * ddSdDesdDes[Si];
        DDIcDDesDDes += psi.dot(dRdS.col(Si)) * ddSdDesdDes[Si];
    }

    return DDIcDDesDDes;
}

MatrixXd adjointDirectHessian(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(n_elem + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(n_elem + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dIcdW
    VectorXd dIcdW(3 * n_elem);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate ddIcdWdW
    SparseMatrix<double> ddIcdWdW;
    ddIcdWdW = evaldIcdWdW(W, dx);
    // Evaluate ddIcdWdDes
    MatrixXd ddIcdWdDes(3 * n_elem, nDesVar);
    ddIcdWdDes = evalddIcdWdS() * dSdDes;
    // Evaluate dIcdS
    VectorXd dIcdS(n_elem + 1);
    dIcdS = evaldIcdS();
    // Evaluate dIcdDes
    VectorXd dIcdDes(nDesVar);
    dIcdDes = dIcdS.transpose() * dSdDes;
    // Evaluate ddIcdDesdDes
    MatrixXd ddIcdDesdDes(nDesVar, nDesVar);
    ddIcdDesdDes = dSdDes.transpose() * evalddIcdSdS() * dSdDes;
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    //// Get Fluxes
    std::vector<double> Flux(3 * (n_elem + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * n_elem, n_elem + 1);
    dRdS = evaldRdS(Flux, area, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * n_elem, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector<double> dt(n_elem, 1);
    SparseMatrix<double> dRdW;
    dRdW = evaldRdW(W, dx, dt, area);
    // Evaluate ddRdWdS
    std::vector <SparseMatrix<double> > ddRdWdS(3 * n_elem);// by (3 * n_elem, n_elem + 1)
    ddRdWdS = evalddRdWdS(W, area);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * n_elem);// by (3 * n_elem, nDesVar)
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
       ddRdWdDes[Ri] = ddRdWdS[Ri] * dSdDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(W, area);

    // *************************************
    // Sparse LU of Jacobian Transpose dRdW.transpose()
    // *************************************
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver1;
    slusolver1.compute(-dRdW.transpose());
    if (slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
    // *************************************
    // Solve for Adjoint (1 Flow Eval)
    // *************************************
    VectorXd psi(3 * n_elem);
    psi = slusolver1.solve(dIcdW);

    // *************************************
    // Evaluate dWdDes (nDesVar Flow Eval)
    // *************************************
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver2;
    slusolver2.compute(-dRdW);
    if (slusolver2.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver2.info()<<std::endl;
    MatrixXd dWdDes(3 * n_elem, nDesVar);
    dWdDes = slusolver2.solve(dRdDes);

    // *************************************
    // Solve for dpsidDes (nDesVar Flow Eval)
    // *************************************
    MatrixXd RHS(3 * n_elem, nDesVar);
    MatrixXd dpsidDes(3 * n_elem, nDesVar);
    RHS = ddIcdWdDes + ddIcdWdW * dWdDes;
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        RHS += psi(Ri) * ddRdWdDes[Ri];
        RHS += psi(Ri) * ddRdWdW[Ri] * dWdDes;
    }
    dpsidDes = slusolver1.solve(RHS);

    // *************************************
    // Evaluate total derivative DDIcDDesDDes
    // *************************************
    MatrixXd DDIcDDesDDes(nDesVar, nDesVar);
    DDIcDDesDDes.setZero();
    DDIcDDesDDes = ddIcdDesdDes
                   + dWdDes.transpose() * ddIcdWdDes
                   + dRdDes.transpose() * dpsidDes;

    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        DDIcDDesDDes += psi(Ri) * (dWdDes.transpose() * ddRdWdDes[Ri]).transpose();
    }
    for (int Si = 0; Si < n_elem + 1; Si++) {
        DDIcDDesDDes += dIcdS(Si) * ddSdDesdDes[Si];
        DDIcDDesDDes += psi.dot(dRdS.col(Si)) * ddSdDesdDes[Si];
    }

    return DDIcDDesDDes;
}
MatrixXd directDirectHessian(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(n_elem + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(n_elem + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dIcdW
    VectorXd dIcdW(3 * n_elem);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate ddIcdWdW
    SparseMatrix<double> ddIcdWdW;
    ddIcdWdW = evaldIcdWdW(W, dx);
    // Evaluate ddIcdWdDes
    MatrixXd ddIcdWdDes(3 * n_elem, nDesVar);
    ddIcdWdDes = evalddIcdWdS() * dSdDes;
    // Evaluate dIcdS
    VectorXd dIcdS(n_elem + 1);
    dIcdS = evaldIcdS();
    // Evaluate dIcdDes
    VectorXd dIcdDes(nDesVar);
    dIcdDes = dIcdS.transpose() * dSdDes;
    // Evaluate ddIcdDesdDes
    MatrixXd ddIcdDesdDes(nDesVar, nDesVar);
    ddIcdDesdDes = dSdDes.transpose() * evalddIcdSdS() * dSdDes;
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    //// Get Fluxes
    std::vector<double> Flux(3 * (n_elem + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * n_elem, n_elem + 1);
    dRdS = evaldRdS(Flux, area, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * n_elem, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector<double> dt(n_elem, 1);
    SparseMatrix<double> dRdW;
    dRdW = evaldRdW(W, dx, dt, area);
    // Evaluate ddRdWdS
    std::vector <SparseMatrix<double> > ddRdWdS(3 * n_elem);// by (3 * n_elem, n_elem + 1)
    ddRdWdS = evalddRdWdS(W, area);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * n_elem);// by (3 * n_elem, nDesVar)
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
       ddRdWdDes[Ri] = ddRdWdS[Ri] * dSdDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(W, area);

    // *************************************
    // Evaluate dWdDes (nDesVar Flow Eval)
    // *************************************
    double errdwddes = -1.0;
    MatrixXd dWdDes(3 * n_elem, nDesVar);
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > factdrdw;
    factdrdw.compute(-dRdW);
    if (factdrdw.info() != 0)
        std::cout<<"Factorization failed. Error: "<<factdrdw.info()<<std::endl;
    if (exactHessian == 1) {
        dWdDes = factdrdw.solve(dRdDes);
    }
    else if (exactHessian == -1 || exactHessian == -3) {
        // Iterative Solution of dWdDes
        dWdDes = solveGMRES(-dRdW,dRdDes);

        std::cout<<"dWdDes ||r|| residual:"<<std::endl;
        std::cout<<(-dRdW*dWdDes - dRdDes).norm()<<std::endl;
        std::cout<<"dWdDes ||r||/||b|| residual:"<<std::endl;
        std::cout<<(-dRdW*dWdDes - dRdDes).norm()/dRdDes.norm()<<std::endl;

        // Direct Solution of dWdDes
        MatrixXd realdWdDes(3*n_elem,nDesVar);
        realdWdDes = factdrdw.solve(dRdDes);

        std::cout<<"Relative error of approximate dWdDes vs exact dWdDes:"<<std::endl;
        std::cout<<(realdWdDes - dWdDes).norm()/realdWdDes.norm()<<std::endl;
        errdwddes = (realdWdDes - dWdDes).norm()/realdWdDes.norm();
    }
    // *************************************
    // Evaluate ddWdDesdDes (nDesVar * (nDesVar+1) / 2 Flow Eval)
    // *************************************
    std::vector <MatrixXd> ddWdDesdDes(3 * n_elem);
    MatrixXd dummy(nDesVar, nDesVar);
    dummy.setZero();
    for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
        ddWdDesdDes[Wi] = dummy;
    }

    if (exactHessian >= -1) {
        VectorXd ddWdDesidDesj(3 * n_elem);
        VectorXd RHS(3 * n_elem);
        for (int di = 0; di < nDesVar; di++) {
            for (int dj = di; dj < nDesVar; dj++) {
                // Evaluate RHS
                RHS.setZero();
                for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
                    RHS[Ri] += dWdDes.col(dj).transpose() * ddRdWdS[Ri] * dSdDes.col(di);
                    RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdS[Ri] * dSdDes.col(dj);
                    RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdW[Ri] * dWdDes.col(dj);
                }

                for (int Si = 0; Si < n_elem + 1; Si++) {
                    RHS += dRdS.col(Si) * ddSdDesdDes[Si](di, dj);
                }

                // Solve
                ddWdDesidDesj = factdrdw.solve(RHS);

                for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
                    ddWdDesdDes[Wi](di, dj) = ddWdDesidDesj[Wi];
                    if (di != dj) {
                        ddWdDesdDes[Wi](dj, di) = ddWdDesdDes[Wi](di, dj);
                    }
                }
            }
        }
    }

    // *************************************
    // Evaluate total derivative DDIcDDesDDes
    // *************************************
    MatrixXd DDIcDDesDDes(nDesVar, nDesVar);
    MatrixXd A(nDesVar, nDesVar);
    A = dWdDes.transpose() * ddIcdWdDes
        + (dWdDes.transpose() * ddIcdWdDes).transpose()
        + dWdDes.transpose() * ddIcdWdW * dWdDes
        + dWdDes.transpose() * ddIcdWdW * dWdDes
        ;
//      + dWdDes.transpose() * ddIcdWdW * dWdDes;
    DDIcDDesDDes = ddIcdDesdDes
                   + dWdDes.transpose() * ddIcdWdDes
                   + (dWdDes.transpose() * ddIcdWdDes).transpose()
                   + dWdDes.transpose() * ddIcdWdW * dWdDes;
    for (int Si = 0; Si < n_elem + 1; Si++) {
        DDIcDDesDDes += dIcdS(Si) * ddSdDesdDes[Si];
    }
    for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
        DDIcDDesDDes += dIcdW(Wi) * ddWdDesdDes[Wi];
    }

    if (exactHessian == -1) {
        for (int i = 0; i<nDesVar; i++) {
            for (int j = 0; j<nDesVar; j++) {
                std::cout << i << "\t" << j 
                    << "\t" << std::setprecision(5) 
                    << errdwddes / (DDIcDDesDDes(i,j)/A(i,j))
                    << std::endl;
            }
        }
    }

    return DDIcDDesDDes;
}

std::vector <MatrixXd> evalddWdDesdDes_FD(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> area,
    std::vector<double> designVar)
{
    std::vector<double> W(3 * n_elem, 0);
    VectorXd W0(3 * n_elem),
             W1(3 * n_elem),
             W2(3 * n_elem),
             W3(3 * n_elem),
             W4(3 * n_elem);

    std::vector <MatrixXd> dwddesddes(3 * n_elem);
    MatrixXd dummy(nDesVar, nDesVar);
    for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
        dwddesddes[Wi] = dummy;
    }

    std::vector<double> tempS(n_elem + 1);

    double dhi, dhj;

    std::vector<double> tempD(nDesVar);

    double h = 1e-3;

    quasiOneD(x, area, W);
    for (int Wi = 0; Wi<3*n_elem; Wi++)
        W0[Wi] = W[Wi];

    for (int i = 0; i < nDesVar; i++) {
        for (int j = i; j < nDesVar; j++) {
            dhi = designVar[i] * h;
            dhj = designVar[j] * h;
            if (i == j)
            {
                tempD = designVar;
                tempD[i] += dhi;
                tempD[j] += dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, tempS, W);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W1[Wi] = W[Wi];
                }

                tempD = designVar;
                tempD[i] += dhi;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, tempS, W);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W2[Wi] = W[Wi];
                }

                tempD = designVar;
                tempD[i] -= dhi;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, tempS, W);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W3[Wi] = W[Wi];
                }

                tempD = designVar;
                tempD[i] -= dhi;
                tempD[j] -= dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, tempS, W);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W4[Wi] = W[Wi];
                }

                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    dwddesddes[Wi](i, j) =
                    (-W1[Wi] + 16*W2[Wi] - 30*W0[Wi] + 16*W3[Wi] - W4[Wi]) / (12 * dhi * dhj);
                }
            }
            else
            {
                tempD = designVar;
                tempD[i] += dhi;
                tempD[j] += dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, tempS, W);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W1[Wi] = W[Wi];
                }

                tempD = designVar;
                tempD[i] += dhi;
                tempD[j] -= dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, tempS, W);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W2[Wi] = W[Wi];
                }

                tempD = designVar;
                tempD[i] -= dhi;
                tempD[j] += dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, tempS, W);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W3[Wi] = W[Wi];
                }

                tempD = designVar;
                tempD[i] -= dhi;
                tempD[j] -= dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, tempS, W);
                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    W4[Wi] = W[Wi];
                }

                for (int Wi = 0; Wi<3*n_elem; Wi++) {
                    dwddesddes[Wi](i, j) = (W1[Wi] - W2[Wi] - W3[Wi] + W4[Wi]) / (4 * dhi * dhj);
                    dwddesddes[Wi](j, i) = dwddesddes[Wi](i, j);
                }

            } // if diag
        }// dj loop
    }// di loop
    return dwddesddes;
}

std::vector <MatrixXd> evalddWdDesdDes(
    std::vector<double> x,
    std::vector<double> dx,
    std::vector<double> W,
    std::vector<double> area,
    std::vector<double> designVar)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(n_elem + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(n_elem + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    // Get Fluxes
    std::vector<double> Flux(3 * (n_elem + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * n_elem, n_elem + 1);
    dRdS = evaldRdS(Flux, area, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * n_elem, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector<double> dt(n_elem, 1);
    SparseMatrix<double> dRdW;
    dRdW = evaldRdW(W, dx, dt, area);
    // Evaluate ddRdWdS
    std::vector <SparseMatrix<double> > ddRdWdS(3 * n_elem);// by (3 * n_elem, n_elem + 1)
    ddRdWdS = evalddRdWdS(W, area);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * n_elem);// by (3 * n_elem, nDesVar)
    for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
        ddRdWdDes[Ri] = ddRdWdS[Ri] * dSdDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix<double> > ddRdWdW(3 * n_elem);// by (3 * n_elem, 3 * n_elem)
    ddRdWdW = evalddRdWdW(W, area);

    // Solve dWdDes
    MatrixXd dWdDes(3 * n_elem, nDesVar);
    // Sparse LU
    SparseLU <SparseMatrix<double>, COLAMDOrdering< int > > slusolver;
    slusolver.analyzePattern(-dRdW);
    slusolver.factorize(-dRdW);
    if (slusolver.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;

    dWdDes = slusolver.solve(dRdDes);

    // Evaluate ddWdDesdDes
    std::vector <MatrixXd> ddWdDesdDes(3 * n_elem);
    MatrixXd dummy(nDesVar, nDesVar);
    VectorXd ddWdDesidDesj(3 * n_elem);
    VectorXd RHS(3 * n_elem);
    for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
        ddWdDesdDes[Wi] = dummy;
    }

    for (int di = 0; di < nDesVar; di++) {
        for (int dj = di; dj < nDesVar; dj++) {
            // Evaluate RHS
            RHS.setZero();
            for (int Ri = 0; Ri < 3 * n_elem; Ri++) {
                RHS[Ri] += dWdDes.col(dj).transpose() * ddRdWdS[Ri] * dSdDes.col(di);
                RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdS[Ri] * dSdDes.col(dj);
                RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdW[Ri] * dWdDes.col(dj);
            }

            for (int Si = 0; Si < n_elem + 1; Si++) {
                RHS += dRdS.col(Si) * ddSdDesdDes[Si](di, dj);
            }

            // Solve
            ddWdDesidDesj = slusolver.solve(RHS);

            for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
                ddWdDesdDes[Wi](di, dj) = ddWdDesidDesj[Wi];
                if (di != dj)
                {
                    ddWdDesdDes[Wi](dj, di) = ddWdDesdDes[Wi](di, dj);
                }
            }
        }
    }

//  std::vector <MatrixXd> ddWdDesdDes_FD = evalddWdDesdDes_FD(x, dx, area, designVar);
//  std::cout<<"Comparing AN and FD ddWdDesdDes"<<std::endl;
//  for (int Wi = 0; Wi < 3 * n_elem; Wi++) {
//      if (std::abs(ddWdDesdDes[Wi].norm()>1e-12))
//      {
//          std::cout<<"Wi "<<Wi<<" diff: "<<
//              (ddWdDesdDes_FD[Wi] - ddWdDesdDes[Wi]).norm()/ddWdDesdDes[Wi].norm()<<std::endl;
//      }
//      else
//      {
//          std::cout<<"Wi "<<Wi<<
//              " zero DD norm. FD norm: "<<(ddWdDesdDes_FD[Wi]).norm()<<std::endl;
//      }
//  }

    return ddWdDesdDes;
}
