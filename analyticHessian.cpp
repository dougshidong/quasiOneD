#include<iostream>
#include<math.h>
#include<vector>
#include<Eigen/Eigen>
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
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar);

std::vector <MatrixXd> evalddWdDesdDes_FD(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar);

MatrixXd directAdjointHessian(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar);
MatrixXd adjointAdjointHessian(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar);
MatrixXd adjointDirectHessian(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar);
MatrixXd directDirectHessian(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar);

MatrixXd getAnalyticHessian(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar,
    int method)
{
    MatrixXd Hessian(nDesVar, nDesVar);

    if(method == 0)
    {
        Hessian = directDirectHessian(x, dx, W, S, designVar);
    }
    else if(method == 1)
    {
        Hessian = adjointDirectHessian(x, dx, W, S, designVar);
    }
    else if(method == 2)
    {
        Hessian = adjointAdjointHessian(x, dx, W, S, designVar);
    }
    else if(method == 3)
    {
        Hessian = directAdjointHessian(x, dx, W, S, designVar);
    }

    return Hessian;
}

MatrixXd directAdjointHessian(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(nx + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dIcdW
    VectorXd dIcdW(3 * nx);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate ddIcdWdW
    SparseMatrix <double> ddIcdWdW;
    ddIcdWdW = evaldIcdWdW(W, dx);
    // Evaluate ddIcdWdDes
    MatrixXd ddIcdWdDes(3 * nx, nDesVar);
    ddIcdWdDes = evalddIcdWdS() * dSdDes;
    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
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
    std::vector <double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    dRdS = evaldRdS(Flux, S, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * nx, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector <double> dt(nx, 1);
    SparseMatrix <double> dRdW;
    dRdW = evaldRdW(W, dx, dt, S);
    // Evaluate ddRdWdS
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);// by (3 * nx, nx + 1)
    ddRdWdS = evalddRdWdS(W, S);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * nx);// by (3 * nx, nDesVar)
    for(int Ri = 0; Ri < 3 * nx; Ri++)
    {
       ddRdWdDes[Ri] = ddRdWdS[Ri] * dSdDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix <double> > ddRdWdW(3 * nx);// by (3 * nx, 3 * nx)
    ddRdWdW = evalddRdWdW(W, S);

    // *************************************
    // Solve for Adjoint (1 Flow Eval)
    // *************************************
    VectorXd psi(3 * nx);
    //SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver1;
    //slusolver1.compute(-dRdW.transpose());
    //if(slusolver1.info() != 0)
    //    std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
    //VectorXd psi(3 * nx);
    //psi = slusolver1.solve(dIcdW);
    psi = solveGMRES(-dRdW.transpose(),dIcdW);

    // *************************************
    // Evaluate dWdDes (nDesVar Flow Eval)
    // *************************************
    MatrixXd dWdDes(3 * nx, nDesVar);
    //SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver2;
    //slusolver2.compute(-dRdW);
    //if(slusolver2.info() != 0)
    //    std::cout<<"Factorization failed. Error: "<<slusolver2.info()<<std::endl;
    //dWdDes = slusolver2.solve(dRdDes);
    dWdDes = solveGMRES(-dRdW,dRdDes);

    // *************************************
    // Evaluate total derivative DDIcDDesDDes
    // *************************************
    MatrixXd DDIcDDesDDes(nDesVar, nDesVar);
    DDIcDDesDDes.setZero();
    DDIcDDesDDes = ddIcdDesdDes;
    DDIcDDesDDes += dWdDes.transpose() * ddIcdWdDes;
    DDIcDDesDDes += (dWdDes.transpose() * ddIcdWdDes).transpose();
    DDIcDDesDDes += dWdDes.transpose() * ddIcdWdW * dWdDes;
    for(int Si = 0; Si < nx + 1; Si++)
    {
        DDIcDDesDDes += dIcdS(Si) * ddSdDesdDes[Si];
        DDIcDDesDDes += psi.dot(dRdS.col(Si)) * ddSdDesdDes[Si];
    }
    for(int Ri = 0; Ri < 3 * nx; Ri++)
    {
//        DDIcDDesDDes += psi(Ri) * ddRdDesdDes;
        DDIcDDesDDes += psi(Ri) * (dWdDes.transpose() * ddRdWdDes[Ri]);
        DDIcDDesDDes += (psi(Ri) * (dWdDes.transpose() * ddRdWdDes[Ri])).transpose();
        DDIcDDesDDes += psi(Ri) * (dWdDes.transpose() * ddRdWdW[Ri] * dWdDes);
    }

    return DDIcDDesDDes;
}
MatrixXd adjointAdjointHessian(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(nx + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dIcdW
    VectorXd dIcdW(3 * nx);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate ddIcdWdW
    SparseMatrix <double> ddIcdWdW;
    ddIcdWdW = evaldIcdWdW(W, dx);
    // Evaluate ddIcdWdDes
    MatrixXd ddIcdWdDes(3 * nx, nDesVar);
    ddIcdWdDes = evalddIcdWdS() * dSdDes;
    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
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
    std::vector <double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    dRdS = evaldRdS(Flux, S, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * nx, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector <double> dt(nx, 1);
    SparseMatrix <double> dRdW;
    dRdW = evaldRdW(W, dx, dt, S);
    // Evaluate ddRdWdS
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);// by (3 * nx, nx + 1)
    ddRdWdS = evalddRdWdS(W, S);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * nx);// by (3 * nx, nDesVar)
    for(int Ri = 0; Ri < 3 * nx; Ri++)
    {
       ddRdWdDes[Ri] = ddRdWdS[Ri] * dSdDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix <double> > ddRdWdW(3 * nx);// by (3 * nx, 3 * nx)
    ddRdWdW = evalddRdWdW(W, S);

    // *************************************
    // Sparse LU of Jacobian Transpose dRdW
    // *************************************
    SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver1;
    slusolver1.compute(-dRdW.transpose());
    if(slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;

    SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver2;
    slusolver2.compute(-dRdW);
    if(slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver2.info()<<std::endl;
    // *************************************
    // Solve for Adjoint 1 psi(1 Flow Eval)
    // *************************************
    VectorXd psi(3 * nx);
    psi = slusolver1.solve(dIcdW);

    // *************************************
    // Solve for Adjoint 2 lambda (nDesVar Flow Eval)
    // *************************************
    MatrixXd lambda(3 * nx, nDesVar);
    lambda = slusolver2.solve(dRdDes);

    // *************************************
    // Solve for Adjoint 3 eta (nDesVar Flow Eval)
    // *************************************
    MatrixXd eta(3 * nx, nDesVar);
    MatrixXd RHS(3 * nx, nDesVar);

    RHS = ddIcdWdDes + ddIcdWdW * lambda;
    for(int Ri = 0; Ri < 3 * nx; Ri++)
    {
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
    for(int Ri = 0; Ri < 3 * nx; Ri++)
    {
//      DDIcDDesDDes += psi(Ri) * ddRdDesdDes[Ri];
        DDIcDDesDDes += psi(Ri) * (lambda.transpose() * ddRdWdDes[Ri]);
    }
    for(int Si = 0; Si < nx + 1; Si++)
    {
        DDIcDDesDDes += dIcdS(Si) * ddSdDesdDes[Si];
        DDIcDDesDDes += psi.dot(dRdS.col(Si)) * ddSdDesdDes[Si];
    }

    return DDIcDDesDDes;
}

MatrixXd adjointDirectHessian(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(nx + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dIcdW
    VectorXd dIcdW(3 * nx);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate ddIcdWdW
    SparseMatrix <double> ddIcdWdW;
    ddIcdWdW = evaldIcdWdW(W, dx);
    // Evaluate ddIcdWdDes
    MatrixXd ddIcdWdDes(3 * nx, nDesVar);
    ddIcdWdDes = evalddIcdWdS() * dSdDes;
    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
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
    std::vector <double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    dRdS = evaldRdS(Flux, S, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * nx, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector <double> dt(nx, 1);
    SparseMatrix <double> dRdW;
    dRdW = evaldRdW(W, dx, dt, S);
    // Evaluate ddRdWdS
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);// by (3 * nx, nx + 1)
    ddRdWdS = evalddRdWdS(W, S);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * nx);// by (3 * nx, nDesVar)
    for(int Ri = 0; Ri < 3 * nx; Ri++)
    {
       ddRdWdDes[Ri] = ddRdWdS[Ri] * dSdDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix <double> > ddRdWdW(3 * nx);// by (3 * nx, 3 * nx)
    ddRdWdW = evalddRdWdW(W, S);

    // *************************************
    // Sparse LU of Jacobian Transpose dRdW.transpose()
    // *************************************
    SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver1;
    slusolver1.compute(-dRdW.transpose());
    if(slusolver1.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver1.info()<<std::endl;
    // *************************************
    // Solve for Adjoint (1 Flow Eval)
    // *************************************
    VectorXd psi(3 * nx);
    psi = slusolver1.solve(dIcdW);

    // *************************************
    // Evaluate dWdDes (nDesVar Flow Eval)
    // *************************************
    SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver2;
    slusolver2.compute(-dRdW);
    if(slusolver2.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver2.info()<<std::endl;
    MatrixXd dWdDes(3 * nx, nDesVar);
    dWdDes = slusolver2.solve(dRdDes);

    // *************************************
    // Solve for dpsidDes (nDesVar Flow Eval)
    // *************************************
    MatrixXd RHS(3 * nx, nDesVar);
    MatrixXd dpsidDes(3 * nx, nDesVar);
    RHS = ddIcdWdDes + ddIcdWdW * dWdDes;
    for(int Ri = 0; Ri < 3 * nx; Ri++)
    {
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
                   + dpsidDes.transpose() * dRdDes;
    for(int Ri = 0; Ri < 3 * nx; Ri++)
    {
        DDIcDDesDDes += psi(Ri) * dWdDes.transpose() * ddRdWdDes[Ri];
    }
    for(int Si = 0; Si < nx + 1; Si++)
    {
        DDIcDDesDDes += dIcdS(Si) * ddSdDesdDes[Si];
        DDIcDDesDDes += psi.dot(dRdS.col(Si)) * ddSdDesdDes[Si];
    }

    return DDIcDDesDDes;
}
MatrixXd directDirectHessian(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(nx + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);

    // *************************************
    // Evaluate Objective Derivatives
    // *************************************
    // Evaluate dIcdW
    VectorXd dIcdW(3 * nx);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate ddIcdWdW
    SparseMatrix <double> ddIcdWdW;
    ddIcdWdW = evaldIcdWdW(W, dx);
    // Evaluate ddIcdWdDes
    MatrixXd ddIcdWdDes(3 * nx, nDesVar);
    ddIcdWdDes = evalddIcdWdS() * dSdDes;
    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
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
    std::vector <double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    dRdS = evaldRdS(Flux, S, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * nx, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector <double> dt(nx, 1);
    SparseMatrix <double> dRdW;
    dRdW = evaldRdW(W, dx, dt, S);
    // Evaluate ddRdWdS
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);// by (3 * nx, nx + 1)
    ddRdWdS = evalddRdWdS(W, S);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * nx);// by (3 * nx, nDesVar)
    for(int Ri = 0; Ri < 3 * nx; Ri++)
    {
       ddRdWdDes[Ri] = ddRdWdS[Ri] * dSdDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix <double> > ddRdWdW(3 * nx);// by (3 * nx, 3 * nx)
    ddRdWdW = evalddRdWdW(W, S);

    // *************************************
    // Evaluate dWdDes (nDesVar Flow Eval)
    // *************************************
    MatrixXd dWdDes(3 * nx, nDesVar);
    SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver;
    slusolver.compute(-dRdW);
    if(slusolver.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;
    dWdDes = slusolver.solve(dRdDes);

    MatrixXd dWdDesGMRES(3 * nx, nDesVar);
    dWdDesGMRES = solveGMRES(-dRdW,dRdDes);

    std::cout<<(dWdDesGMRES-dWdDes).norm()/dWdDes.norm()<<std::endl;
    // *************************************
    // Evaluate ddWdDesdDes (nDesVar * (nDesVar+1) / 2 Flow Eval)
    // *************************************
    std::vector <MatrixXd> ddWdDesdDes(3 * nx);
    MatrixXd dummy(nDesVar, nDesVar);
    for(int Wi = 0; Wi < 3 * nx; Wi++)
    {
        ddWdDesdDes[Wi] = dummy;
    }

    VectorXd ddWdDesidDesj(3 * nx);
    VectorXd RHS(3 * nx);
    for(int di = 0; di < nDesVar; di++)
    {
        for(int dj = di; dj < nDesVar; dj++)
        {
            // Evaluate RHS
            RHS.setZero();
            for(int Ri = 0; Ri < 3 * nx; Ri++)
            {
                RHS[Ri] += dWdDes.col(dj).transpose() * ddRdWdS[Ri] * dSdDes.col(di);
                RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdS[Ri] * dSdDes.col(dj);
                RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdW[Ri] * dWdDes.col(dj);
            }

            for(int Si = 0; Si < nx + 1; Si++)
            {
                RHS += dRdS.col(Si) * ddSdDesdDes[Si](di, dj);
            }

            // Solve
            ddWdDesidDesj = slusolver.solve(RHS);

            for(int Wi = 0; Wi < 3 * nx; Wi++)
            {
                ddWdDesdDes[Wi](di, dj) = ddWdDesidDesj[Wi];
                if(di != dj)
                {
                    ddWdDesdDes[Wi](dj, di) = ddWdDesdDes[Wi](di, dj);
                }
            }
        }
    }

    // *************************************
    // Evaluate total derivative DDIcDDesDDes
    // *************************************
    MatrixXd DDIcDDesDDes(nDesVar, nDesVar);
    DDIcDDesDDes = ddIcdDesdDes
                   + dWdDes.transpose() * ddIcdWdDes
                   + (dWdDes.transpose() * ddIcdWdDes).transpose()
                   + dWdDes.transpose() * ddIcdWdW * dWdDes;
    for(int Si = 0; Si < nx + 1; Si++)
    {
        DDIcDDesDDes += dIcdS(Si) * ddSdDesdDes[Si];
    }
    for(int Wi = 0; Wi < 3 * nx; Wi++)
    {
        DDIcDDesDDes += dIcdW(Wi) * ddWdDesdDes[Wi];
    }

    return DDIcDDesDDes;
}

std::vector <MatrixXd> evalddWdDesdDes_FD(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar)
{
    std::vector <double> W(3 * nx, 0);
    VectorXd W0(3 * nx),
             W1(3 * nx),
             W2(3 * nx),
             W3(3 * nx),
             W4(3 * nx);

    std::vector <MatrixXd> dwddesddes(3 * nx);
    MatrixXd dummy(nDesVar, nDesVar);
    for(int Wi = 0; Wi < 3 * nx; Wi++)
    {
        dwddesddes[Wi] = dummy;
    }

    std::vector <double> tempS(nx + 1);

    double dhi, dhj;

    std::vector <double> tempD(nDesVar);

    double h = 1e-3;

    quasiOneD(x, dx, S, W);
    for(int Wi = 0; Wi<3*nx; Wi++)
        W0[Wi] = W[Wi];

    for(int i = 0; i < nDesVar; i++)
    {
        for(int j = i; j < nDesVar; j++)
        {
            dhi = designVar[i] * h;
            dhj = designVar[j] * h;
            if(i == j)
            {
                tempD = designVar;
                tempD[i] += dhi;
                tempD[j] += dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W1[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] += dhi;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W2[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] -= dhi;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W3[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] -= dhi;
                tempD[j] -= dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W4[Wi] = W[Wi];

                for(int Wi = 0; Wi<3*nx; Wi++)
                    dwddesddes[Wi](i, j) =
                    (-W1[Wi] + 16*W2[Wi] - 30*W0[Wi] + 16*W3[Wi] - W4[Wi]) / (12 * dhi * dhj);
            }
            else
            {
                tempD = designVar;
                tempD[i] += dhi;
                tempD[j] += dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W1[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] += dhi;
                tempD[j] -= dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W2[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] -= dhi;
                tempD[j] += dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W3[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] -= dhi;
                tempD[j] -= dhj;
                tempS = evalS(tempD, x, dx, desParam);
                quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W4[Wi] = W[Wi];

                for(int Wi = 0; Wi<3*nx; Wi++)
                {
                    dwddesddes[Wi](i, j) = (W1[Wi] - W2[Wi] - W3[Wi] + W4[Wi]) / (4 * dhi * dhj);
                    dwddesddes[Wi](j, i) = dwddesddes[Wi](i, j);
                }

            } // if diag
        }// dj loop
    }// di loop
    return dwddesddes;
}

std::vector <MatrixXd> evalddWdDesdDes(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar)
{
    // *************************************
    // Evaluate Area to Design Derivatives
    // *************************************
    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(nx + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);
    // *************************************
    // Evaluate Residual Derivatives
    // *************************************
    // Get Fluxes
    std::vector <double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);
    // Evaluate dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    dRdS = evaldRdS(Flux, S, W);
    // Evaluate dRdDes
    MatrixXd dRdDes(3 * nx, nDesVar);
    dRdDes = dRdS * dSdDes;
    // Evaluate dRdW
    std::vector <double> dt(nx, 1);
    SparseMatrix <double> dRdW;
    dRdW = evaldRdW(W, dx, dt, S);
    // Evaluate ddRdWdS
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);// by (3 * nx, nx + 1)
    ddRdWdS = evalddRdWdS(W, S);
    // Evaluate ddRdWdDes
    std::vector <MatrixXd> ddRdWdDes(3 * nx);// by (3 * nx, nDesVar)
    for(int Ri = 0; Ri < 3 * nx; Ri++)
    {
        ddRdWdDes[Ri] = ddRdWdS[Ri] * dSdDes;
    }
    // Evaluate ddRdWdW
    std::vector < SparseMatrix <double> > ddRdWdW(3 * nx);// by (3 * nx, 3 * nx)
    ddRdWdW = evalddRdWdW(W, S);

    // Solve dWdDes
    MatrixXd dWdDes(3 * nx, nDesVar);
    // Sparse LU
    SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver;
    slusolver.analyzePattern(-dRdW);
    slusolver.factorize(-dRdW);
    if(slusolver.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;

    dWdDes = slusolver.solve(dRdDes);

    // Evaluate ddWdDesdDes
    std::vector <MatrixXd> ddWdDesdDes(3 * nx);
    MatrixXd dummy(nDesVar, nDesVar);
    VectorXd ddWdDesidDesj(3 * nx);
    VectorXd RHS(3 * nx);
    for(int Wi = 0; Wi < 3 * nx; Wi++)
    {
        ddWdDesdDes[Wi] = dummy;
    }

    for(int di = 0; di < nDesVar; di++)
    {
        for(int dj = di; dj < nDesVar; dj++)
        {
            // Evaluate RHS
            RHS.setZero();
            for(int Ri = 0; Ri < 3 * nx; Ri++)
            {
                RHS[Ri] += dWdDes.col(dj).transpose() * ddRdWdS[Ri] * dSdDes.col(di);
                RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdS[Ri] * dSdDes.col(dj);
                RHS[Ri] += dWdDes.col(di).transpose() * ddRdWdW[Ri] * dWdDes.col(dj);
            }

            for(int Si = 0; Si < nx + 1; Si++)
            {
                RHS += dRdS.col(Si) * ddSdDesdDes[Si](di, dj);
            }

            // Solve
            ddWdDesidDesj = slusolver.solve(RHS);

            for(int Wi = 0; Wi < 3 * nx; Wi++)
            {
                ddWdDesdDes[Wi](di, dj) = ddWdDesidDesj[Wi];
                if(di != dj)
                {
                    ddWdDesdDes[Wi](dj, di) = ddWdDesdDes[Wi](di, dj);
                }
            }
        }
    }

//  std::vector <MatrixXd> ddWdDesdDes_FD = evalddWdDesdDes_FD(x, dx, S, designVar);
//  std::cout<<"Comparing AN and FD ddWdDesdDes"<<std::endl;
//  for(int Wi = 0; Wi < 3 * nx; Wi++)
//  {
//      if(std::abs(ddWdDesdDes[Wi].norm()>1e-12))
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
