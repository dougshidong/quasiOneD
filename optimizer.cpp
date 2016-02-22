#include<iostream>
#include<math.h>
#include"quasiOneD.h"
#include"grid.h"
#include"adjoint.h"
#include"globals.h"
#include<vector>
#include<iomanip>
#include<Eigen/Dense>
#include <complex>
#include<stdlib.h>//exit
std::vector <std::complex<double> > finiteD(std::vector <std::complex<double> > x,
                             std::vector <std::complex<double> > dx,
                             std::vector <std::complex<double> > S, 
                             std::vector <std::complex<double> > designVar, 
                             std::complex<double> h, 
                             std::complex<double> currentI);

std::complex<double> stepBacktrackUncons(std::vector <std::complex<double> > designVar, 
                           std::vector <std::complex<double> > pk, 
                           std::vector <std::complex<double> > gradient, 
                           std::complex<double> currentI, 
                           std::vector <std::complex<double> > x, 
                           std::vector <std::complex<double> > dx);

std::vector <std::complex<double> > BFGS(std::vector <std::complex<double> > oldH, 
                          std::vector <std::complex<double> > gradList, 
                          std::vector <std::complex<double> > searchD);

void design(std::vector <std::complex<double> > x, std::vector <std::complex<double> > dx,
            std::vector <std::complex<double> > S, std::vector <std::complex<double> > designVar)
{
    std::vector <std::complex<double> > W(3 * nx, 0);

    
    std::vector <std::complex<double> > normGradList, gradList;
    std::vector <std::complex<double> > H(nDesVar * nDesVar, 0);
    std::complex<double> normGrad = 1;
    std::complex<double> tolGrad = 1e-10;
    std::complex<double> currentI;

    int maxDesign = 10000;

    int printConv = 1;

    std::complex<double> alpha = 1;

    std::complex<double> h = 1e-25;
    std::vector <std::complex<double> > gradient(nDesVar), gradientFD(nDesVar),
                pk(nDesVar),
                searchD(nDesVar),
                dgradient(nDesVar);

    int rc;
    std::vector <std::complex<double> > psi(3 * nx, 0);
    
    quasiOneD(x, dx, S, designVar, W);
    
    currentI = quasiOneD(x, dx, S, designVar, W);
    gradient = adjoint(x, dx, S, W, psi, designVar);
    gradientFD = finiteD(x, dx, S, designVar, h, currentI);

    std::cout<<"Gradient Relative Error:"<<std::endl;
    for(int i = 0; i < nDesVar; i++)
    {
        double gradientDiff = fabs(std::real((gradient[i] - gradientFD[i]) / (gradient[i])));
        std::cout<<gradientDiff<<std::endl;
    }

    
//    gradient = gradientFD;
      exit(EXIT_FAILURE); 
    // Initialize B
    for(int r = 0; r < nDesVar; r++)
    for(int c = 0; c < nDesVar; c++)
    {
        if(r == c)
        {
            rc = r * nDesVar + c;
            H[rc] = 1;
        }
        
    }

    normGradList.push_back(1);
    int iDesign = 0;

    // Design Loop
    while(std::real(normGrad) > std::real(tolGrad) && iDesign < maxDesign)
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
        S = evalS(designVar, x, dx);
        currentI = quasiOneD(x, dx, S, designVar, W);

        gradientFD = finiteD(x, dx, S, designVar, h, currentI);
        gradient = adjoint(x, dx, S, W, psi, designVar);
//        gradient = gradientFD;

        for(int i = 0; i < nDesVar; i++)
        {
            gradList.push_back(gradient[i]);
        }

//      1  =  Steepest Descent
//      3  =  Newton
//      4  =  Quasi-Newton (BFGS)
        if(descentType == 1)
        {
            for(int i = 0; i < nDesVar; i++)
                pk[i] =  - gradient[i];
        }
        if(descentType == 4)
        {
            if(iDesign > 1)
            {
                H = BFGS(H, gradList, searchD);
            }

            for(int r = 0; r < nDesVar; r++)
            {
                pk[r] = 0;
                for(int c = 0; c<nDesVar; c++)
                {
                    rc = r * nDesVar + c;
                    pk[r] +=  - H[rc] * gradient[c];
                }
            }
        }
        std::cout<<"pk:\n";
        for(int i = 0; i < nDesVar; i++)
            std::cout<<pk[i]<<std::endl;

                alpha = stepBacktrackUncons(designVar, pk, gradient, currentI,
                      x, dx);
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
        
        std::cout<<"End of Design Iteration: "<<iDesign<<std::endl<<std::endl<<std::endl;
    }

    std::cout<<"Final Gradient:"<<std::endl;
    for(int i = 0; i < nDesVar; i++)
        std::cout<<gradient[i]<<std::endl;

    std::cout<<std::endl<<"Final Design:"<<std::endl;
    for(int i = 0; i < nDesVar; i++)
        std::cout<<designVar[i]<<std::endl;


    S = evalS(designVar, x, dx);

    std::cout<<"Fitness: "<<quasiOneD(x, dx, S, designVar, W)<<std::endl;





    return;



}

std::vector <std::complex<double> > BFGS(std::vector <std::complex<double> > oldH,
                          std::vector <std::complex<double> > gradList,
                          std::vector <std::complex<double> > searchD)
{
    int ls = gradList.size();
    int rc;
    std::vector <std::complex<double> > newH(nDesVar * nDesVar);
    Eigen::VectorXcd dg(nDesVar), dx(nDesVar);
    Eigen::MatrixXcd dH(nDesVar, nDesVar), a(nDesVar, nDesVar), b(nDesVar, nDesVar);

    for(int i = 0; i < nDesVar; i++)
    {
        dg(i) = gradList[ls - nDesVar + i] - gradList[ls - 2 * nDesVar + i];
        dx(i) = searchD[i];
    }
    
    Eigen::Map <Eigen::Matrix<std::complex<double>, - 1, - 1, Eigen::RowMajor> > 
        cH(oldH.data(), nDesVar, nDesVar);

    a = ((dx.transpose() * dg + dg.transpose() * cH * dg)(0) * (dx * dx.transpose()))
         / ((dx.transpose() * dg)(0) * (dx.transpose() * dg)(0));
    b = (cH * dg * dx.transpose() + dx * dg.transpose() * cH) / (dx.transpose() * dg)(0);

    dH = a - b;

    std::cout<<"Current Inverse Hessian:"<<std::endl;
    for(int r = 0; r < nDesVar; r++)
    {
        std::cout<<"\n";
        for(int c = 0; c < nDesVar; c++)
        {
            rc = r * nDesVar + c;
            newH[rc] = oldH[rc] + dH(r, c);
            std::cout<<newH[rc]<<"\t\t";
        }
    }

    std::cout<<"\n\n";

    return newH;
}

std::vector <std::complex<double> > finiteD(std::vector <std::complex<double> > x, 
                             std::vector <std::complex<double> > dx,
                             std::vector <std::complex<double> > S, 
                             std::vector <std::complex<double> > designVar, 
                             std::complex<double> h, 
                             std::complex<double> currentI)
{
    std::vector <std::complex<double> > W(3 * nx, 0);

    // Method
    // 1  =  Forward
    // 2  =  Backward
    // 3  =  Central
    std::vector <std::complex<double> > grad(nx);
    std::vector <std::complex<double> > tempS(nx + 1);
    
    std::complex<double> I0, I1, I2, dh;

    std::vector <std::complex<double> > tempD(nDesVar);

    if(std::real(currentI) < 0 && gradientType != 3)
    {
        I0 = quasiOneD(x, dx, S, designVar, W);

        //std::cout<<"I0 = "<<std::setprecision(15)<<I0<<std::endl;
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
            
            //for(int k = 0; k < nDesVar; k++)
            //    std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;

            tempS = evalS(tempD, x, dx);

            I1 = quasiOneD(x, dx, tempS, tempD, W);

            grad[i] = (I1 - I0) / dh;
//          std::cout<<"I1 = "<<std::setprecision(15)<<I1<<std::endl;
    
        }
        else if(gradientType == 2)
        {
            tempD[i] -= dh;
        
            //for(int k = 0; k < nDesVar; k++)
            //    std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;

            tempS = evalS(tempD, x, dx);

            I2 = quasiOneD(x, dx, tempS, tempD, W);

            grad[i] = (I0 - I2) / dh;
            //std::cout<<"I2 = "<<std::setprecision(15)<<I2<<std::endl;

        }
        else if(gradientType == 3)
        {
            tempD[i] += dh;
            //for(int k = 0; k < nDesVar; k++)
            //    std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;


            tempS = evalS(tempD, x, dx);


            I1 = quasiOneD(x, dx, tempS, tempD, W);

            //std::cout<<"I1 = "<<std::setprecision(15)<<I1<<std::endl;

            tempD = designVar;
            tempD[i] -= dh;

            //for(int k = 0; k < nDesVar; k++)
            //    std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;

            tempS = evalS(tempD, x, dx);

            I2 = quasiOneD(x, dx, tempS, tempD, W);

            //std::cout<<"I2 = "<<std::setprecision(15)<<I2<<std::endl;

            grad[i] = (I1 - I2) / (2.0 * dh);
        }
        else if(gradientType == 4)
        {
            tempD[i] += dh * std::complex<double>(0.0,1.0);

            tempS = evalS(tempD, x, dx);

            I1 = quasiOneD(x, dx, tempS, tempD, W);

            std::cout<<"doing complex step"<<std::endl;

            grad[i] = std::imag(I1) / (dh);
        }
    }
    std::cout<<"Gradient from FD: "<<std::endl;
    for(int i = 0; i < 3; i++)
        std::cout<<grad[i]<<std::endl;

    return grad;
}



std::complex<double> stepBacktrackUncons(std::vector <std::complex<double> > designVar,
                           std::vector <std::complex<double> > pk,
                           std::vector <std::complex<double> > gradient, 
                           std::complex<double> currentI,
                           std::vector <std::complex<double> > x, 
                           std::vector <std::complex<double> > dx)
{
    std::vector <std::complex<double> > W(3 * nx, 0);

    std::complex<double> alpha = 1;
    std::complex<double> c1 = 1e-4;
    std::vector <std::complex<double> > tempS(nx + 1);
    std::complex<double> newVal;

    std::complex<double> c_pk_grad = 0;

    std::vector <std::complex<double> > tempD(nDesVar);

    for(int i = 0; i < nDesVar; i++)
    {
        c_pk_grad += gradient[i] * pk[i];
    }

    c_pk_grad *= c1;
    
    for(int i = 0; i < nDesVar; i++)
    {
        tempD[i] = designVar[i] + alpha * pk[i];
    }

    tempS = evalS(tempD, x, dx);
    newVal = quasiOneD(x, dx, tempS, tempD, W);


    while(std::real(newVal) > std::real((currentI + alpha * c_pk_grad)) && std::real(alpha) > 0.0001)
    {
        alpha = alpha * 0.5;
        std::cout<<"Alpha Reduction: "<<alpha<<std::endl;
    
        for(int i = 0; i < nDesVar; i++)
            tempD[i] = designVar[i] + alpha * pk[i];
        tempS = evalS(tempD, x, dx);
        newVal = quasiOneD(x, dx, tempS, tempD, W);

    }

    return alpha;
}   

