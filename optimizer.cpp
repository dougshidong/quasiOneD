#include<iostream>
#include<fstream>
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
#include<stdlib.h>     /* srand, rand */

using namespace Eigen;

MatrixXd finiteD2g(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar,
    double h);
MatrixXd finiteD2(std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar,
    double h,
    double currentI,
    int &possemidef);

double stepBacktrackUncons(
    double alpha,
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

VectorXd implicitSmoothing(VectorXd gradient, double epsilon);

void test_grad(int gradientType1, int gradientType2, double currentI, 
	std::vector<double> x, std::vector<double> dx, 
	std::vector<double> S, std::vector<double> W, 
	std::vector<double> designVar, VectorXd psi)
{
	printf("Comparing GradientType1: %d and GradientType2: %d", gradientType1, gradientType2);
	int n_des = designVar.size();
    VectorXd gradient1(n_des);
    VectorXd gradient2(n_des);
    gradient1 = getGradient(gradientType1, currentI, x, dx, S, W, designVar, psi);
    gradient2 = getGradient(gradientType2, currentI, x, dx, S, W, designVar, psi);

	printf("GradientType1: %d \t GradientType2: %d \t Relative Difference\n", gradientType1, gradientType2);
	for (int i = 0; i < n_des; i++) {
		double g1 = gradient1[i];
		double g2 = gradient2[i];
		printf("%23.15e \t %23.15e \t %23.15e\n", g1, g2, (g1-g2)/(g1+1e-15));
	}
	return;
}
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
    std::vector <double> svdvalues;
    std::vector <double> svdvaluesreal;
    std::vector <double> Hcond;
    MatrixXd H(nDesVar, nDesVar), H_BFGS(nDesVar, nDesVar), realH(nDesVar, nDesVar);
    double normGrad;
    double currentI;
    double alpha;

    int printConv = 1;

    VectorXd pk(nDesVar), searchD(nDesVar);

    clock_t tic = clock();
    clock_t toc;
    double elapsed;

    std::ofstream myfile;
    myfile.open("convergence.dat");
    myfile << " Iteration \t Cost Function \t Gradient Norm \t Average Error \n";

    quasiOneD(x, dx, S, W);
    currentI = evalFitness(dx, W);

    VectorXd psi(3 * nx);

    VectorXd gradient(nDesVar);
    VectorXd oldGrad(nDesVar); //BFGS
    gradient = getGradient(gradientType, currentI, x, dx, S, W, designVar, psi);

	bool testingGradient = true;
	testingGradient = false;
	if (testingGradient) {
		test_grad(gradientType, -3, currentI, x, dx, S, W, designVar, psi);
		exit(0);
	}


//  exactHessian = 1; // Calculate exact Hessian to compare with BFGS
//  realH = getAnalyticHessian(x, dx, W, S, designVar, 3);
//  int tempi;
//  double err;

//  H= finiteD2g(x, dx, S, designVar, 1.0e-07);
//  err = (realH - H).norm()/realH.norm();
//  std::cout<<"Exact - FDg: "<<err<<std::endl;

//  H= getAnalyticHessian(x, dx, W, S, designVar, 1);
//  err = (realH - H).norm()/realH.norm();
//  std::cout<<"DA - AD: "<<err<<std::endl;
//  std::cout<<std::setprecision(15)<<H<<std::endl;

//  H= getAnalyticHessian(x, dx, W, S, designVar, 2);
//  err = (realH - H).norm()/realH.norm();
//  std::cout<<"DA - AA: "<<err<<std::endl;
//  std::cout<<std::setprecision(15)<<H<<std::endl;


//  H= getAnalyticHessian(x, dx, W, S, designVar, 0);
//  err = (realH - H).norm()/realH.norm();
//  std::cout<<"DA - DD: "<<err<<std::endl;
//  std::cout<<std::setprecision(15)<<H<<std::endl;
//  return;

    // Initialize B
    H.setIdentity();
    H = H * 1.0;
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

//          std::cout<<"Current Shape:\n";
//          for(int i = 0; i < nx + 1; i++)
//              std::cout<<S[i]<<std::endl;
        }
        std::cout<<"Current Fitness: "<<currentI<<std::endl;

//      1  =  Steepest Descent
//      2  =  Quasi-Newton (BFGS)
//      3  =  Newton
//      4  =  Truncated Newton with Adjoint-Direct Matrix-Vector Product
        if(descentType == 1)
        {
            //int expo = rand() % 5 + 1 - 3;
            //double averageerr = 0;
            //for(int i = 0; i<nDesVar; i++){
            //    srand (time(NULL));
            //    double fMin = -2.0;
            //    double fMax = 2.0;
            //    double expo = (double)rand() / RAND_MAX;
            //    expo =  fMin + expo * (fMax - fMin);
            //    expo =  0.0;
            //    pk(i) =  -10*gradient(i)*pow(10,expo);
            //    averageerr += fabs(expo)/nDesVar;
            //    //pk(i) =  -gradient(i);
            //}
            //myfile << iDesign << "\t" << currentI <<"\t"<< normGrad << "\t" << averageerr << "\n";
            //myfile.flush();
            pk =  -500*gradient;
        }
        else if(descentType == 2)
        {
            if(iDesign > 1)
            {
                H_BFGS = BFGS(H, oldGrad, gradient, searchD);
                H = H_BFGS;
            }

//          realH = getAnalyticHessian(x, dx, W, S, designVar, 2);
//          JacobiSVD<MatrixXd> svd1(H.inverse(), ComputeFullU | ComputeFullV);
//          JacobiSVD<MatrixXd> svd2(realH, ComputeFullU | ComputeFullV);

//          std::cout<<"svd1"<<std::endl;
//          std::cout<<svd1.singularValues()<<std::endl;
//          std::cout<<"svd2"<<std::endl;
//          std::cout<<svd2.singularValues()<<std::endl;
//          for(int i = 0; i < nDesVar; i++)
//          {
//              svdvalues.push_back(svd1.singularValues()(i));
//              svdvaluesreal.push_back(svd2.singularValues()(i));
//          }
//
//          std::cout<<"svd singular values error"<<std::endl;
//          std::cout<<
//              (svd1.singularValues()-svd2.singularValues()).norm()
//              /svd2.singularValues().norm()<<std::endl;

//          std::cout<<"svd singular vectors error"<<std::endl;
//          std::cout<<
//              (svd1.matrixV()-svd2.matrixV()).norm()
//              /svd2.matrixV().norm()<<std::endl;
            // Eigenvalues are not returned in ascending order.
//          std::cout<<"eig1"<<std::endl;
//          std::cout<<H.inverse().eigenvalues()<<std::endl;
//          std::cout<<"eig2"<<std::endl;
//          std::cout<<realH.eigenvalues()<<std::endl;
//          std::cout<<"Eig error"<<std::endl;
//          std::cout<<
//              (H.inverse().eigenvalues()-realH.eigenvalues()).norm()
//              /realH.eigenvalues().norm()<<std::endl;

//          realH = realH.inverse();
//          Hcond.push_back(checkCond(realH.inverse()));
//          double err = (realH - H).norm()/realH.norm();
//          std::cout<<"Hessian error: "<<err<<std::endl;
//          Herror.push_back(err);

            pk = -H * gradient;
        }
        else if(descentType == 3)
        {
            H = getAnalyticHessian(x, dx, W, S, designVar, hessianType);
            realH = getAnalyticHessian(x, dx, W, S, designVar, 2);
            double err = (realH - H).norm()/realH.norm();
            std::cout<<"Hessian error: "<<err<<std::endl;
            Herror.push_back(err);

            checkCond(H);
            H = invertHessian(H);
            Hcond.push_back(checkCond(realH.inverse()));

            pk = -H * gradient;
        }
//        std::cout<<realH<<std::endl;
        
//      std::cout<<"pk before smoothing:\n"<<std::endl;
//      std::cout<<pk<<std::endl;

//      pk = implicitSmoothing(pk, 0.5);
        alpha = 1.0;

        std::cout<<"gradient:\n"<<std::endl;
        std::cout<<-gradient<<std::endl;
        std::cout<<"pk:\n"<<std::endl;
        std::cout<<pk<<std::endl;

        currentI = stepBacktrackUncons(alpha, designVar, searchD, pk, gradient, currentI, x, dx, W);

        S = evalS(designVar, x, dx, desParam);
        oldGrad = gradient;
        gradient = getGradient(gradientType, currentI, x, dx, S, W, designVar, psi);

        normGrad = 0;
        for(int i = 0; i < nDesVar; i++)
            normGrad += pow(gradient[i], 2);
        normGrad = sqrt(normGrad);
        normGradList.push_back(normGrad);

        toc = clock();
        elapsed = (double)(toc-tic) / CLOCKS_PER_SEC;
        timeVec.push_back(elapsed);
        std::cout<<"Time: "<<elapsed<<std::endl;

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
    outVec("svd.dat", "w", svdvalues);
    outVec("svdreal.dat", "w", svdvaluesreal);

    myfile.close();

    return;
}

double stepBacktrackUncons(
    double alpha,
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

    double c1 = 1e-4;
    std::vector <double> tempS(nx + 1);
    double newVal;

    double c_pk_grad = 0;


    c_pk_grad = c1 * gradient.dot(pk);

    std::vector <double> tempD(nDesVar);
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

MatrixXd finiteD2g(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> S,
    std::vector <double> designVar,
    double h)
{
    MatrixXd Hessian(nDesVar, nDesVar);
    std::vector <double> W(3 * nx, 0);
    std::vector <double> tempS(nx + 1);
    std::vector <double> tempD(nDesVar);
    VectorXd gradp(nDesVar);
    VectorXd gradn(nDesVar);
    VectorXd psi(3 * nx);

    double currentI = -1;
    Hessian.setZero();
    for(int i = 0; i < nDesVar; i++)
    {
        tempD = designVar;
        tempD[i] += h;
        tempS = evalS(tempD, x, dx, desParam);
        quasiOneD(x, dx, tempS, W);
        gradp = getGradient(gradientType, currentI, x, dx, tempS, W, tempD, psi);

        tempD = designVar;
        tempD[i] -= h;
        tempS = evalS(tempD, x, dx, desParam);
        quasiOneD(x, dx, tempS, W);
        gradn = getGradient(gradientType, currentI, x, dx, tempS, W, tempD, psi);
        for(int j = 0; j < nDesVar; j++)
        {
            Hessian(i, j) += (gradp(j) - gradn(j)) / (2*h);
            Hessian(j, i) += (gradp(j) - gradn(j)) / (2*h);
        }
    }
    Hessian = Hessian / 2.0;
    return Hessian;
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
    Hessian = Hessian + Hessian.transpose();
    Hessian = Hessian / 2.0;
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

    dg = currentg - oldg;
    dx = searchD;

    a = ((dx.transpose() * dg + dg.transpose() * oldH * dg)(0) * (dx * dx.transpose()))
         / ((dx.transpose() * dg)(0) * (dx.transpose() * dg)(0));
    b = (oldH * dg * dx.transpose() + dx * dg.transpose() * oldH) / (dx.transpose() * dg)(0);

    dH = a - b;

    newH = oldH + dH;

    return newH;
}

VectorXd implicitSmoothing(VectorXd gradient, double epsilon)
{
    int n = gradient.size();
    MatrixXd A(n, n);
    A.setZero();

    for(int i = 0; i < n - 1; i++)
    {
        A(i  , i) = 1.0 + 2.0 * epsilon;
        A(i+1, i) = -epsilon;
        A(i, i+1) = -epsilon;
    }
    A(n-1,n-1) = 1.0 + 2.0 * epsilon;

    LLT<MatrixXd> llt;
    llt.compute(A);
    if(llt.info() != 0)
        std::cout<<"Factorization failed. Error: "<<llt.info()<<std::endl;
    VectorXd smoothGrad = llt.solve(gradient);

    return smoothGrad;
}
