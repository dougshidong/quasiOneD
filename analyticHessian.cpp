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

using namespace Eigen;

std::vector <MatrixXd> evalddWdDesdDes(
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

    Hessian = directDirectHessian(x, dx, W, S, designVar);

    return Hessian;
}

MatrixXd directDirectHessian(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> W,
    std::vector <double> S,
    std::vector <double> designVar)
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

    // Evaluate dIcdW
    VectorXd dIcdW(3 * nx);
    dIcdW = evaldIcdW(W, dx);
    // Evaluate ddIcdWdW
    SparseMatrix <double> ddIcdWdW;
    ddIcdWdW = evaldIcdWdW(W, dx);

    // Evaluate dIcdS
    VectorXd dIcdS(nx + 1);
    dIcdS.setZero();
    // Evaluate ddIcdSdS
    MatrixXd ddIcdSdS(nx + 1, nx + 1);
    ddIcdSdS = evalddIcdSdS();

    // Evaluate ddIcdWdS
    MatrixXd ddIcdWdS(3 * nx, nx + 1);
    ddIcdWdS = evalddIcdWdS();

    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, designVar.size());
    dSdDes = evaldSdDes(x, dx, designVar);
    // Evaluate ddSdDesdDes
    std::vector <MatrixXd> ddSdDesdDes(nx + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);

    // Evaluate dWdDes
    MatrixXd dWdDes(3 * nx, nDesVar);
    dWdDes = evaldWdDes(x, dx, S, W, designVar);
    // Evaluate ddWdDesdDes
    std::vector <MatrixXd> ddWdDesdDes(3 * nx);
    ddWdDesdDes = evalddWdDesdDes(x, dx, W, S, designVar);

    MatrixXd ddIcdDesdDes(nDesVar, nDesVar);
    ddIcdDesdDes.setZero();
    double a, b, c, d;
    for(int di = 0; di < nDesVar; di++)
    {
        for(int dj = di; dj < nDesVar; dj++)
        {
            a = dSdDes.col(di).transpose() * ddIcdSdS * dSdDes.col(dj); // Zero
            b = dWdDes.col(dj).transpose() * ddIcdWdS * dSdDes.col(di); // Zero
            c = dWdDes.col(di).transpose() * ddIcdWdS * dSdDes.col(dj); // Zero
            d = dWdDes.col(di).transpose() * ddIcdWdW * dWdDes.col(dj);
            ddIcdDesdDes(di, dj) = a + b + c + d;

//          d = dWdDes.col(di).transpose() * ddIcdWdW * dWdDes.col(dj);
//          ddIcdDesdDes(di, dj) = d;

            for(int Si = 0; Si < nx + 1; Si++)
            {
                  ddIcdDesdDes(di, dj) += dIcdS(Si) * ddSdDesdDes[Si](di, dj); // Zero
            }
            for(int Wi = 0; Wi < 3 * nx; Wi++)
            {
                ddIcdDesdDes(di, dj) += dIcdW(Wi) * ddWdDesdDes[Wi](di, dj);
            }

            if(di != dj)
            {
                ddIcdDesdDes(dj, di) = ddIcdDesdDes(di, dj);
            }
        }
    }

    return ddIcdDesdDes;
}

std::vector <MatrixXd> evaldwddesddes(
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

    double I, I1, I2, I3, I4, dhi, dhj;

    std::vector <double> tempD(nDesVar);

    double h = 1e-4;

    I1 = quasiOneD(x, dx, S, W);
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
                I1 = quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W1[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] += dhi;
                tempS = evalS(tempD, x, dx, desParam);
                I2 = quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W2[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] -= dhi;
                tempS = evalS(tempD, x, dx, desParam);
                I3 = quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W3[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] -= dhi;
                tempD[j] -= dhj;
                tempS = evalS(tempD, x, dx, desParam);
                I4 = quasiOneD(x, dx, tempS, W);
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
                I1 = quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W1[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] += dhi;
                tempD[j] -= dhj;
                tempS = evalS(tempD, x, dx, desParam);
                I2 = quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W2[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] -= dhi;
                tempD[j] += dhj;
                tempS = evalS(tempD, x, dx, desParam);
                I3 = quasiOneD(x, dx, tempS, W);
                for(int Wi = 0; Wi<3*nx; Wi++)
                    W3[Wi] = W[Wi];

                tempD = designVar;
                tempD[i] -= dhi;
                tempD[j] -= dhj;
                tempS = evalS(tempD, x, dx, desParam);
                I4 = quasiOneD(x, dx, tempS, W);
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
    //Get Primitive Variables
    std::vector <double> rho(nx), u(nx), e(nx);
    std::vector <double> T(nx), p(nx), c(nx), Mach(nx);
    WtoP(W, rho, u, e, p, c, T);

    // ****************
    // First Derivative
    // ****************

    // Get Fluxes
    std::vector <double> Flux(3 * (nx + 1), 0);
    getFlux(Flux, W);

    // Evaluate dRdS
    MatrixXd dRdS(3 * nx, nx + 1);
    dRdS = evaldRdS(Flux, S, W);

    // Evaluate dRdW
    std::vector <double> dt(nx, 1);
    SparseMatrix <double> dRdW;
    dRdW = evaldRdW(W, dx, dt, S, u[0]/c[0]);

    // Evaluate dSdDes
    MatrixXd dSdDes(nx + 1, nDesVar);
    dSdDes = evaldSdDes(x, dx, designVar);

    //Evaluate dRdDes
    MatrixXd dRdDes(3 * nx, nDesVar);
    dRdDes = dRdS * dSdDes;

    // Solve dWdDes
    MatrixXd dWdDes(3 * nx, nDesVar);
    // Sparse LU
    SparseLU <SparseMatrix <double>, COLAMDOrdering< int > > slusolver;
    slusolver.analyzePattern(-dRdW);
    slusolver.factorize(-dRdW);
    if(slusolver.info() != 0)
        std::cout<<"Factorization failed. Error: "<<slusolver.info()<<std::endl;

    dWdDes = slusolver.solve(dRdDes);

    // ****************
    // Second Derivative
    // ****************

    // Evaluate ddRdWdS // Confirmed with FD
    std::vector <SparseMatrix <double> > ddRdWdS(3 * nx);// by (3 * nx, nx + 1)
    ddRdWdS = evalddRdWdS(W, S);

    // Evaluate ddRdWdW // Confirmed with FD
    std::vector < SparseMatrix <double> > ddRdWdW(3 * nx);// by (3 * nx, 3 * nx)
    ddRdWdW = evalddRdWdW(W, S);

    // Evaluate ddSdDesdDes // Confirmed with FD
    std::vector <MatrixXd> ddSdDesdDes(nx + 1); // by (nDes, nDes)
    ddSdDesdDes = evalddSdDesdDes(x, dx, designVar);

    // RHS
    VectorXd RHS(3 * nx);
    // ddWdDesidDesj
    VectorXd ddWdDesidDesj(3 * nx);

    // Evaluate ddWdDesdDes
    std::vector <MatrixXd> ddWdDesdDes(3 * nx);
    MatrixXd dummy(nDesVar, nDesVar);
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

    std::vector <MatrixXd> ddWdDesdDes_FD = evaldwddesddes(x, dx, S, designVar);

//  std::cout<<"Comparing AN and FD ddWdDesdDes"<<std::endl;
//  for(int Wi = 0; Wi < 3 * nx; Wi++)
//  {
//      std::cout<<(ddWdDesdDes_FD[Wi] - ddWdDesdDes[Wi]).norm()/ddWdDesdDes[Wi].norm()<<std::endl;
//  }

    return ddWdDesdDes;
}
