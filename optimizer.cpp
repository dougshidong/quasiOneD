#include<iostream>
#include<math.h>
#include<vector>
#include<iomanip>
#include<Eigen/Dense>
#include<stdlib.h>//exit
#include"quasiOneD.h"
#include"fitness.h"
#include"grid.h"
#include"adjoint.h"
#include"directDifferentiation.h"
#include"globals.h"
#include"gradient.h"
#include"analyticHessian.h"
#include"output.h"
#include<time.h>

using namespace Eigen;

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

void checkCond(MatrixXd H);
MatrixXd invertHessian(MatrixXd H);
LLT<MatrixXd> checkPosDef(MatrixXd H);
LLT<MatrixXd> makePosDef(MatrixXd H);

void design(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar)
{
    std::vector <double> W(3 * nx, 0);

    std::vector <double> normGradList, gradList;
    std::vector <double> timeVec;
    std::vector <double> Herror;
    MatrixXd H(nDesVar, nDesVar), H_BFGS(nDesVar, nDesVar), realH(nDesVar, nDesVar);
    double normGrad;
    double currentI;

    int maxDesign = 10000;

    int printConv = 1;

    double alpha = 1;

    VectorXd pk(nDesVar), searchD(nDesVar);

    clock_t tic = clock();
    clock_t toc;
    double elapsed;


    quasiOneD(x, dx, S, W);
    currentI = evalFitness(dx, W);

    VectorXd gradient(nDesVar);
    gradient = getGradient(gradientType, currentI, x, dx, S, W, designVar);

    // Initialize B
    H.setIdentity();
    if(exactHessian == 1)
    {
        H = getAnalyticHessian(x, dx, W, S, designVar, hessianType);
//      H = finiteD2(x, dx, S, designVar, h, currentI, possemidef);
        checkCond(H);
        H = invertHessian(H);
    }

    normGrad = 0;
    for(int i = 0; i < nDesVar; i++)
        normGrad += pow(gradient[i], 2);
    normGrad = sqrt(normGrad);
    int iDesign = 0;

    // Design Loop
    while(normGrad > gradConv && iDesign < maxDesign)
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
        quasiOneD(x, dx, S, W);
        currentI = evalFitness(dx, W);
        std::cout<<"Current Fitness: "<<currentI<<std::endl;

        gradient = getGradient(gradientType, currentI, x, dx, S, W, designVar);

        for(int i = 0; i < nDesVar; i++)
        {
            gradList.push_back(gradient[i]);
        }

//      1  =  Steepest Descent
//      2  =  Quasi-Newton (BFGS)
//      3  =  Newton
        if(descentType == 1)
        {
            for(int i = 0; i < nDesVar; i++)
                pk[i] =  - gradient[i];
        }
        else if(descentType == 2)
        {
            if(iDesign > 1)
            {
                H_BFGS = BFGS(H, gradList, searchD);
                H = H_BFGS;
            }

            realH = getAnalyticHessian(x, dx, W, S, designVar, hessianType).inverse();
            Herror.push_back((realH - H).norm()/realH.norm());

//          std::cout<<"Current Inverse Hessian:"<<std::endl;
//          std::cout<<H<<std::endl;
//          std::cout<<"\n\n";

            pk = -H * gradient;
        }
        else if(descentType == 3)
        {
            if(iDesign > 1)
            {
                H = getAnalyticHessian(x, dx, W, S, designVar, hessianType);
                checkCond(H);
                H = invertHessian(H);
            }

            std::cout<<"Current Inverse Hessian:"<<std::endl;
            std::cout<<H<<std::endl;
            std::cout<<"\n\n";

            pk = -H * gradient;
        }
        std::cout<<"pk:\n";
        for(int i = 0; i < nDesVar; i++)
            std::cout<<pk[i]<<std::endl;

        alpha = stepBacktrackUncons(designVar, pk, gradient, currentI, x, dx);
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
            normGrad += pow(gradient[i], 2);
        normGrad = sqrt(normGrad);
        normGradList.push_back(normGrad);

        toc = clock();
        elapsed = (double)(toc-tic) / CLOCKS_PER_SEC;
        timeVec.push_back(elapsed);

        std::cout<<"End of Design Iteration: "<<iDesign<<std::endl<<std::endl<<std::endl;
    }

    std::cout<<"Final Gradient:"<<std::endl;
    for(int i = 0; i < nDesVar; i++)
        std::cout<<gradient[i]<<std::endl;

    std::cout<<std::endl<<"Final Design:"<<std::endl;
    for(int i = 0; i < nDesVar; i++)
        std::cout<<designVar[i]<<std::endl;


    S = evalS(designVar, x, dx, desParam);

    quasiOneD(x, dx, S, W);
    std::cout<<"Fitness: "<<evalFitness(dx, W)<<std::endl;

    outVec("OptConv.dat", "w", normGradList);
    outVec("OptTime.dat", "w", timeVec);
    outVec("HessianErr.dat", "w", Herror);


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

    return newH;
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
    quasiOneD(x, dx, tempS, W);
    newVal = evalFitness(dx, W);


    while(newVal > (currentI + alpha/2.0 * c_pk_grad) && alpha > 1e-16)
    {
        alpha = alpha * 0.5;
        std::cout<<"Alpha Reduction: "<<alpha<<std::endl;

        for(int i = 0; i < nDesVar; i++)
            tempD[i] = designVar[i] + alpha * pk[i];
        tempS = evalS(tempD, x, dx, desParam);
        quasiOneD(x, dx, tempS, W);
        newVal = evalFitness(dx, W);
        std::cout<<"newVal: "<<newVal<<std::endl;
        std::cout<<"currentI + alpha/2.0 * c_pk_grad: "<<
        currentI + alpha/ 2.0 * c_pk_grad<<std::endl;
    }
    if(alpha < 1e-16) std::cout<<"Error. Can't find step size"<<std::endl;

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
        quasiOneD(x, dx, S, W);
        I = evalFitness(dx, W);
    }
    else
    {
        I = currentI;
    }
    for(int i = 0; i < nDesVar; i++)
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
            I1 = evalFitness(dx, W);

            tempD = designVar;
            tempD[i] += dhi;
            tempS = evalS(tempD, x, dx, desParam);
            quasiOneD(x, dx, tempS, W);
            I2 = evalFitness(dx, W);

            tempD = designVar;
            tempD[i] -= dhi;
            tempS = evalS(tempD, x, dx, desParam);
            quasiOneD(x, dx, tempS, W);
            I3 = evalFitness(dx, W);

            tempD = designVar;
            tempD[i] -= dhi;
            tempD[j] -= dhj;
            tempS = evalS(tempD, x, dx, desParam);
            quasiOneD(x, dx, tempS, W);
            I4 = evalFitness(dx, W);
            Hessian(i, j) = (-I1 + 16*I2 - 30*I + 16*I3 - I4) / (12 * dhi * dhj);
        }
        else
        {
            tempD = designVar;
            tempD[i] += dhi;
            tempD[j] += dhj;
            tempS = evalS(tempD, x, dx, desParam);
            quasiOneD(x, dx, tempS, W);
            I1 = evalFitness(dx, W);

            tempD = designVar;
            tempD[i] += dhi;
            tempD[j] -= dhj;
            tempS = evalS(tempD, x, dx, desParam);
            quasiOneD(x, dx, tempS, W);
            I2 = evalFitness(dx, W);

            tempD = designVar;
            tempD[i] -= dhi;
            tempD[j] += dhj;
            tempS = evalS(tempD, x, dx, desParam);
            quasiOneD(x, dx, tempS, W);
            I3 = evalFitness(dx, W);

            tempD = designVar;
            tempD[i] -= dhi;
            tempD[j] -= dhj;
            tempS = evalS(tempD, x, dx, desParam);
            quasiOneD(x, dx, tempS, W);
            I4 = evalFitness(dx, W);

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

void checkCond(MatrixXd H)
{
    JacobiSVD<MatrixXd> svd(H);
    double svdmax = svd.singularValues()(0);
    double svdmin = svd.singularValues()(svd.singularValues().size()-1);
    double cond = svdmax / svdmin;
    std::cout<<"Condition Number of H:"<<std::endl;
    std::cout<<cond<<std::endl;
}

MatrixXd invertHessian(MatrixXd H)
{
    LLT<MatrixXd> llt = checkPosDef(H);
    return llt.solve(MatrixXd::Identity(H.rows(), H.rows()));
}

LLT<MatrixXd> checkPosDef(MatrixXd H)
{
    LLT<MatrixXd> llt;
    VectorXcd eigval = H.eigenvalues();
    double shift = 1e-5;
    if(eigval.real().minCoeff() < 0)
    {
        MatrixXd eye(H.rows(),H.rows());
        eye.setIdentity();
        std::cout<<"Matrix is not Positive Semi-Definite"<<std::endl;
        std::cout<<"Eigenvalues:"<<std::endl;
        std::cout<<eigval<<std::endl;
        //llt = makePosDef(H);
        llt.compute(H + (shift - eigval.real().minCoeff()) * eye);
        checkCond(H + (shift - eigval.real().minCoeff()) * eye);
    }
    else
    {
        llt.compute(H);
    }
    return llt;
}

LLT<MatrixXd> makePosDef(MatrixXd H)
{
    int n = H.rows();
    MatrixXd modH(n,n);
    MatrixXd eye(n,n);
    eye.setIdentity();

    LLT<MatrixXd> llt;

    double beta = H.norm();
    double tau;

    if(H.diagonal().minCoeff() > 0.0) tau = 0.0;
    else tau = beta / 2.0;

    llt.compute(H + tau * eye);
    while(llt.info() != Success)
    {
        std::cout<<"Tried adjusting by tau = "<<tau<<std::endl;
        tau = std::max(2.0 * tau, beta / 2.0);
        llt.compute(H + tau * eye);
    }
    std::cout<<"Adjusted by tau = "<<tau<<std::endl;

    return llt;
}
