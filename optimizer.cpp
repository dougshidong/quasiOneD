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
    std::vector <double> &designVar,
    VectorXd &searchD,
    VectorXd pk,
    VectorXd gradient,
    double currentI,
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> &W);

MatrixXd BFGS(
    MatrixXd oldH,
    VectorXd oldg,
    VectorXd currentg,
    VectorXd searchD);

double checkCond(MatrixXd H);
MatrixXd invertHessian(MatrixXd H);
LLT<MatrixXd> checkPosDef(MatrixXd H);

void design(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar)
{
    std::vector <double> W(3 * nx, 0);

    std::vector <double> normGradList;
    std::vector <double> timeVec;
    std::vector <double> Herror;
    std::vector <double> Hcond;
    MatrixXd H(nDesVar, nDesVar), H_BFGS(nDesVar, nDesVar), realH(nDesVar, nDesVar);
    double normGrad;
    double currentI;

    int printConv = 1;

    VectorXd pk(nDesVar), searchD(nDesVar);

    clock_t tic = clock();
    clock_t toc;
    double elapsed;


    quasiOneD(x, dx, S, W);
    currentI = evalFitness(dx, W);

    VectorXd gradient(nDesVar);
    VectorXd oldGrad(nDesVar); //BFGS
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
        std::cout<<"Current Fitness: "<<currentI<<std::endl;

//      1  =  Steepest Descent
//      2  =  Quasi-Newton (BFGS)
//      3  =  Newton
        if(descentType == 1) pk =  -gradient;
        else if(descentType == 2)
        {
            if(iDesign > 1)
            {
                H_BFGS = BFGS(H, oldGrad, gradient, searchD);
                H = H_BFGS;
            }

            realH = getAnalyticHessian(x, dx, W, S, designVar, hessianType).inverse();
            Hcond.push_back(checkCond(realH.inverse()));
            double err = (realH - H).norm()/realH.norm();
            std::cout<<"Hessian error: "<<err<<std::endl;
            Herror.push_back(err);

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
            realH = getAnalyticHessian(x, dx, W, S, designVar, 2).inverse();
            Hcond.push_back(checkCond(realH.inverse()));
            double err = (realH - H).norm()/realH.norm();
            std::cout<<"Hessian error: "<<err<<std::endl;
            Herror.push_back(err);

            pk = -H * gradient;
        }
//      std::cout<<"Current Inverse Hessian:"<<std::endl;
//      std::cout<<H<<std::endl;
//      std::cout<<"\n\n";

        std::cout<<"pk:\n"<<std::endl;
        std::cout<<pk<<std::endl;

        currentI = stepBacktrackUncons(designVar, searchD, pk, gradient, currentI, x, dx, W);

        S = evalS(designVar, x, dx, desParam);
        oldGrad = gradient;
        gradient = getGradient(gradientType, currentI, x, dx, S, W, designVar);

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
    std::cout<<gradient<<std::endl;

    std::cout<<std::endl<<"Final Design:"<<std::endl;
    for(int i = 0; i < nDesVar; i++)
        std::cout<<designVar[i]<<std::endl;

    std::cout<<"Fitness: "<<evalFitness(dx, W)<<std::endl;

    outVec("OptConv.dat", "w", normGradList);
    outVec("OptTime.dat", "w", timeVec);
    outVec("HessianErr.dat", "w", Herror);
    outVec("HessianCond.dat", "w", Hcond);


    return;
}

double stepBacktrackUncons(
    std::vector <double> &designVar,
    VectorXd &searchD,
    VectorXd pk,
    VectorXd gradient,
    double currentI,
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> &W)
{
//  std::vector <double> W(3 * nx, 0);

    double alpha = 1;
    double c1 = 1e-4;
    std::vector <double> tempS(nx + 1);
    double newVal;

    double c_pk_grad = 0;

    std::vector <double> tempD(nDesVar);

    c_pk_grad = c1 * gradient.dot(pk);

    for(int i = 0; i < nDesVar; i++)
    {
        tempD[i] = designVar[i] + alpha * pk[i];
    }

    tempS = evalS(tempD, x, dx, desParam);
    quasiOneD(x, dx, tempS, W);
    newVal = evalFitness(dx, W);

    while(newVal > (currentI + alpha * c_pk_grad) && alpha > 1e-16)
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

    designVar = tempD;
    for(int i = 0; i < nDesVar; i++)
    {
        searchD[i] = alpha * pk[i];
    }

    return newVal;
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
    return Hessian;
}

double checkCond(MatrixXd H)
{
    JacobiSVD<MatrixXd> svd(H);
    double svdmax = svd.singularValues()(0);
    double svdmin = svd.singularValues()(svd.singularValues().size()-1);
    double cond = svdmax / svdmin;
    std::cout<<"Condition Number of H:"<<std::endl;
    std::cout<<cond<<std::endl;

    return cond;
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
        llt.compute(H + (shift - eigval.real().minCoeff()) * eye);
        checkCond(H + (shift - eigval.real().minCoeff()) * eye);
    }
    else
    {
        llt.compute(H);
    }
    return llt;
}

MatrixXd BFGS(
    MatrixXd oldH,
    VectorXd oldg,
    VectorXd currentg,
    VectorXd searchD)
{
    MatrixXd newH(nDesVar, nDesVar);
    VectorXd dg(nDesVar), dx(nDesVar);
    MatrixXd dH(nDesVar, nDesVar), a(nDesVar, nDesVar), b(nDesVar, nDesVar);

    dg = oldg - currentg;
    dx = searchD;

    a = ((dx.transpose() * dg + dg.transpose() * oldH * dg)(0) * (dx * dx.transpose()))
         / ((dx.transpose() * dg)(0) * (dx.transpose() * dg)(0));
    b = (oldH * dg * dx.transpose() + dx * dg.transpose() * oldH) / (dx.transpose() * dg)(0);

    dH = a - b;

    newH = oldH + dH;

    return newH;
}
