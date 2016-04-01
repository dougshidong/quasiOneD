#include<iostream>
#include<math.h>
#include<vector>
#include<iomanip>
#include<Eigen/Dense>
#include<stdlib.h>//exit
#include"quasiOneD.h"
#include"grid.h"
#include"adjoint.h"
#include"directDifferentiation.h"
#include"globals.h"
#include"analyticHessian.h"

using namespace Eigen;
VectorXd finiteD(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar,
    double h,
    double currentI);

MatrixXd finiteD2(std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar,
    double h,
    double currentI,
    int &possemidef);

double stepBacktrackUncons(
    std::vector <double> designVar,
    VectorXd pk,
    VectorXd gradient,
    double currentI,
    std::vector <double> x,
    std::vector <double> dx);

MatrixXd BFGS(
    MatrixXd oldH,
    std::vector <double> gradList,
    VectorXd searchD);

void design(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar)
{
    std::vector <double> W(3 * nx, 0);


    std::vector <double> normGradList, gradList;
    MatrixXd H(nDesVar, nDesVar), H_FD(nDesVar, nDesVar), H_BFGS(nDesVar, nDesVar);
    double normGrad = 1;
    double tolGrad = 1e-10;
    double currentI;
    int possemidef;

    int maxDesign = 10000;

    int printConv = 1;

    double alpha = 1;

    double h = 1e-8;

    VectorXd pk(nDesVar), searchD(nDesVar), dgradient(nDesVar);

    std::vector <double> psi(3 * nx, 0);

    quasiOneD(x, dx, S, W);

    currentI = quasiOneD(x, dx, S, W);

    VectorXd gradientAV(nDesVar), gradientFD(nDesVar), gradientDD(nDesVar);
    gradientAV = adjoint(x, dx, S, W, psi, designVar);
    gradientDD = directDifferentiation(x, dx, S, W, designVar);
    gradientFD = finiteD(x, dx, S, designVar, h, currentI);

//  std::cout<<"Gradient Relative Error:"<<std::endl;
//  std::cout<<"Norm(AV - DD)/Norm(AV): ";
//  std::cout<<(gradientAV - gradientDD).norm() / gradientAV.norm()<<std::endl;
//  std::cout<<"Norm(AV - FD)/Norm(AV): ";
//  std::cout<<(gradientAV - gradientFD).norm() / gradientAV.norm()<<std::endl;
//  std::cout<<"Norm(DD - FD)/Norm(DD): ";
//  std::cout<<(gradientDD - gradientFD).norm() / gradientDD.norm()<<std::endl;

    int Hmethod = 0;

    // Initialize B
    for(int r = 0; r < nDesVar; r++)
    for(int c = 0; c < nDesVar; c++)
    {
        if(r == c)
        {
            H(r, c) = 1;
        }
    }
    MatrixXd HDD(nDesVar, nDesVar), 
             HAD(nDesVar, nDesVar),
             HAA(nDesVar, nDesVar),
             HDA(nDesVar, nDesVar),
             HFD(nDesVar, nDesVar);
    
    HDD = getAnalyticHessian(x, dx, W, S, designVar, 0);
    std::cout<<"DD Hessian"<<std::endl;
    std::cout<<HDD<<std::endl;
    
    HAD = getAnalyticHessian(x, dx, W, S, designVar, 1);
    std::cout<<"AD Hessian"<<std::endl;
    std::cout<<HAD<<std::endl;
    
    HAA = getAnalyticHessian(x, dx, W, S, designVar, 2);
    std::cout<<"AA Hessian"<<std::endl;
    std::cout<<HAA<<std::endl;

    HDA = getAnalyticHessian(x, dx, W, S, designVar, 3);
    std::cout<<"DA Hessian"<<std::endl;
    std::cout<<HDA<<std::endl;

    HFD = finiteD2(x, dx, S, designVar, h, currentI, possemidef);
    std::cout<<"FD Hessian"<<std::endl;
    std::cout<<HFD<<std::endl;

    std::cout<<"Norm(DD - AD)/Norm(DD): ";
    std::cout<<(HDD - HAD).norm() / HDD.norm()<<std::endl;
    std::cout<<"Norm(DD - AA)/Norm(DD): ";
    std::cout<<(HDD - HAA).norm() / HDD.norm()<<std::endl;
    std::cout<<"Norm(DD - DA)/Norm(DD): ";
    std::cout<<(HDD - HDA).norm() / HDD.norm()<<std::endl;
    std::cout<<"Norm(DD - FD)/Norm(DD): ";
    std::cout<<(HDD - HFD).norm() / HDD.norm()<<std::endl;
    exit(EXIT_FAILURE);

    normGradList.push_back(1);
    int iDesign = 0;

    // Design Loop
    while(normGrad > tolGrad && iDesign < maxDesign)
    {
        iDesign++ ;

        if(printConv == 1)
        {
            std::cout<<"Iteration :"<<iDesign<<
                "    GradientNorm: "<<normGrad<<std::endl;
            std::cout<<"Current Design:\n";
            for(int i = 0; i < nDesVar; i++)
                std::cout<<designVar[i]<<std::endl;


        }
        S = evalS(designVar, x, dx, desParam);
        currentI = quasiOneD(x, dx, S, W);

        gradientFD = finiteD(x, dx, S, designVar, h, currentI);
        gradientAV = adjoint(x, dx, S, W, psi, designVar);

        for(int i = 0; i < nDesVar; i++)
        {
            gradList.push_back(gradientAV[i]);
        }

//      1  =  Steepest Descent
//      3  =  Newton
//      4  =  Quasi-Newton (BFGS)
        if(descentType == 1)
        {
            for(int i = 0; i < nDesVar; i++)
                pk[i] =  - gradientAV[i];
        }
        if(descentType == 4)
        {
            if(iDesign > 1)
            {
                H_BFGS = BFGS(H, gradList, searchD);
                H_FD = finiteD2(x, dx, S, designVar, h, currentI, possemidef).inverse();
                H = H_BFGS;
//              if(possemidef == 1) H = H_FD;
//              H = H_FD;
            }

            for(int r = 0; r < nDesVar; r++)
            {
                pk[r] = 0;
                for(int c = 0; c<nDesVar; c++)
                {
                    pk[r] +=  - H(r, c) * gradientAV[c];
                }
            }
        }
        std::cout<<"pk:\n";
        for(int i = 0; i < nDesVar; i++)
            std::cout<<pk[i]<<std::endl;

        alpha = stepBacktrackUncons(designVar, pk, gradientAV, currentI, x, dx);
        std::cout<<"Alpha = "<<alpha<<std::endl;
        for(int i = 0; i < nDesVar; i++)
        {
            searchD[i] = alpha * pk[i];
        }

        std::cout<<"Search Direction: :"<<std::endl;
        for(int i = 0; i < nDesVar; i++)
            std::cout<<searchD[i]<<std::endl;


        for(int i = 0; i < nDesVar; i++)
            designVar[i] = designVar[i] + searchD[i];

        normGrad = 0;
        for(int i = 0; i < nDesVar; i++)
            normGrad += pow(gradientAV[i], 2);
        normGrad = sqrt(normGrad);
        normGradList.push_back(normGrad);

        std::cout<<"End of Design Iteration: "<<iDesign<<std::endl<<std::endl<<std::endl;
    }

    std::cout<<"Final Gradient:"<<std::endl;
    for(int i = 0; i < nDesVar; i++)
        std::cout<<gradientAV[i]<<std::endl;

    std::cout<<std::endl<<"Final Design:"<<std::endl;
    for(int i = 0; i < nDesVar; i++)
        std::cout<<designVar[i]<<std::endl;


    S = evalS(designVar, x, dx, desParam);

    std::cout<<"Fitness: "<<quasiOneD(x, dx, S, W)<<std::endl;


    return;
}

MatrixXd BFGS(
    MatrixXd oldH,
    std::vector <double> gradList,
    VectorXd searchD)
{
    int ls = gradList.size();
    Eigen::MatrixXd newH(nDesVar, nDesVar);
    Eigen::VectorXd dg(nDesVar), dx(nDesVar);
    Eigen::MatrixXd dH(nDesVar, nDesVar), a(nDesVar, nDesVar), b(nDesVar, nDesVar);

    for(int i = 0; i < nDesVar; i++)
    {
        dg(i) = gradList[ls - nDesVar + i] - gradList[ls - 2 * nDesVar + i];
        dx(i) = searchD[i];
    }

    a = ((dx.transpose() * dg + dg.transpose() * oldH * dg)(0) * (dx * dx.transpose()))
         / ((dx.transpose() * dg)(0) * (dx.transpose() * dg)(0));
    b = (oldH * dg * dx.transpose() + dx * dg.transpose() * oldH) / (dx.transpose() * dg)(0);

    dH = a - b;

    newH = oldH + dH;
    std::cout<<"Current Inverse Hessian:"<<std::endl;
    std::cout<<newH<<std::endl;
    std::cout<<"\n\n";

    return newH;
}

VectorXd finiteD(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar,
    double h,
    double currentI)
{
    std::vector <double> W(3 * nx, 0);

    // Method
    // 1  =  Forward
    // 2  =  Backward
    // 3  =  Central
    VectorXd grad(nDesVar);
    std::vector <double> tempS(nx + 1);

    double I0, I1, I2, dh;

    std::vector <double> tempD(nDesVar);

    if(currentI < 0 && gradientType != 3)
    {
        I0 = quasiOneD(x, dx, S, W);
    }
    else
    {
        I0 = currentI;
    }
    for(int i = 0; i < nDesVar; i++)
    {

        dh = designVar[i] * h;
        tempD = designVar;

        if(gradientType == 1)
        {
            tempD[i] += dh;
            tempS = evalS(tempD, x, dx, desParam);
            I1 = quasiOneD(x, dx, tempS, W);
            grad[i] = (I1 - I0) / dh;
        }
        else if(gradientType == 2)
        {
            tempD[i] -= dh;
            tempS = evalS(tempD, x, dx, desParam);
            I2 = quasiOneD(x, dx, tempS, W);
            grad[i] = (I0 - I2) / dh;
        }
        else if(gradientType == 3)
        {
            tempD[i] += dh;
            tempS = evalS(tempD, x, dx, desParam);
            I1 = quasiOneD(x, dx, tempS, W);
            tempD = designVar;

            tempD[i] -= dh;
            tempS = evalS(tempD, x, dx, desParam);
            I2 = quasiOneD(x, dx, tempS, W);
            grad[i] = (I1 - I2) / (2 * dh);
        }
    }
    std::cout<<"Gradient from FD: "<<std::endl;
    for(int i = 0; i < nDesVar; i++)
        std::cout<<grad[i]<<std::endl;

    return grad;
}

double stepBacktrackUncons(
    std::vector <double> designVar,
    VectorXd pk,
    VectorXd gradient,
    double currentI,
    std::vector <double> x,
    std::vector <double> dx)
{
    std::vector <double> W(3 * nx, 0);

    double alpha = 1;
    double c1 = 1e-4;
    std::vector <double> tempS(nx + 1);
    double newVal;

    double c_pk_grad = 0;

    std::vector <double> tempD(nDesVar);

    for(int i = 0; i < nDesVar; i++)
    {
        c_pk_grad += pow(gradient[i] * pk[i], 2);
    }

//    c_pk_grad *= c1;
    c_pk_grad = sqrt(c_pk_grad);

    for(int i = 0; i < nDesVar; i++)
    {
        tempD[i] = designVar[i] + alpha * pk[i];
    }

    tempS = evalS(tempD, x, dx, desParam);
    newVal = quasiOneD(x, dx, tempS, W);


    while(newVal > (currentI + alpha/2.0 * c_pk_grad) && alpha > 1e-20)
    {
        alpha = alpha * 0.5;
        std::cout<<"Alpha Reduction: "<<alpha<<std::endl;

        for(int i = 0; i < nDesVar; i++)
            tempD[i] = designVar[i] + alpha * pk[i];
        tempS = evalS(tempD, x, dx, desParam);
        newVal = quasiOneD(x, dx, tempS, W);
        std::cout<<"newVal: "<<newVal<<std::endl;
        std::cout<<"currentI + alpha/2.0 * c_pk_grad: "<<
        currentI + alpha/ 2.0 * c_pk_grad<<std::endl;
    }

    return alpha;
}

MatrixXd finiteD2(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar,
    double h,
    double currentI,
    int &possemidef)
{
    std::vector <double> W(3 * nx, 0);
    MatrixXd Hessian(nDesVar, nDesVar);
    std::vector <double> tempS(nx + 1);

    double I, I1, I2, I3, I4, dhi, dhj;

    std::vector <double> tempD(nDesVar);

    h = 1e-4;

    if(currentI < 0 && gradientType != 3)
    {
        I = quasiOneD(x, dx, S, W);
    }
    else
    {
        I = currentI;
    }
    for(int i = 0; i < nDesVar; i++)
    for(int j = i; j < nDesVar; j++)
    {
        std::cout<<"i = "<<i<<", j = "<<j<<"  ";
        dhi = designVar[i] * h;
        dhj = designVar[j] * h;
        if(i == j)
        {
            tempD = designVar;
            tempD[i] += dhi;
            tempD[j] += dhj;
            tempS = evalS(tempD, x, dx, desParam);
            I1 = quasiOneD(x, dx, tempS, W);

            tempD = designVar;
            tempD[i] += dhi;
            tempS = evalS(tempD, x, dx, desParam);
            I2 = quasiOneD(x, dx, tempS, W);

            tempD = designVar;
            tempD[i] -= dhi;
            tempS = evalS(tempD, x, dx, desParam);
            I3 = quasiOneD(x, dx, tempS, W);

            tempD = designVar;
            tempD[i] -= dhi;
            tempD[j] -= dhj;
            tempS = evalS(tempD, x, dx, desParam);
            I4 = quasiOneD(x, dx, tempS, W);
            Hessian(i, j) = (-I1 + 16*I2 - 30*I + 16*I3 - I4) / (12 * dhi * dhj);
        }
        else
        {
            tempD = designVar;
            tempD[i] += dhi;
            tempD[j] += dhj;
            tempS = evalS(tempD, x, dx, desParam);
            I1 = quasiOneD(x, dx, tempS, W);

            tempD = designVar;
            tempD[i] += dhi;
            tempD[j] -= dhj;
            tempS = evalS(tempD, x, dx, desParam);
            I2 = quasiOneD(x, dx, tempS, W);

            tempD = designVar;
            tempD[i] -= dhi;
            tempD[j] += dhj;
            tempS = evalS(tempD, x, dx, desParam);
            I3 = quasiOneD(x, dx, tempS, W);

            tempD = designVar;
            tempD[i] -= dhi;
            tempD[j] -= dhj;
            tempS = evalS(tempD, x, dx, desParam);
            I4 = quasiOneD(x, dx, tempS, W);

            Hessian(i, j) = (I1 - I2 - I3 + I4) / (4 * dhi * dhj);
            Hessian(j, i) = Hessian(i, j);
        }
    }
//  VectorXcd eigval = Hessian.eigenvalues();
//  if(eigval.real().minCoeff() < 0)
//  {
//      MatrixXd eye = MatrixXd(Hessian.rows(), Hessian.cols()).setIdentity();
//      std::cout<<"Matrix is not Positive Semi-Definite"<<std::endl;
//      std::cout<<Hessian<<std::endl;
//      std::cout<<eigval<<std::endl;
//      Hessian = Hessian + fabs(eigval.real().minCoeff()) * eye
//                        + 0.001 * eye;
//      std::cout<<"New Eigenvalues:"<<Hessian.eigenvalues()<<std::endl;
//      possemidef = 0;
//  }
//  else
//  {
//      possemidef = 1;
//  }

//  std::cout<<"Inverse Hessian from FD: "<<std::endl;
//  std::cout<<Hessian.inverse()<<std::endl;


//  JacobiSVD<MatrixXd> svd(Hessian);
//  double svdmax = svd.singularValues()(0);
//  double svdmin = svd.singularValues()(svd.singularValues().size()-1);
//  double cond = svdmax / svdmin;
//  std::cout<<"Condition Number Hessian"<<std::endl;
//  std::cout<<cond<<std::endl;
    return Hessian;
}
