// Calculates the discrete costate fluxes
#include<iostream>
#include<math.h>
#include"quasiOneD.h"
#include<vector>
#include<Eigen/Dense>
#include <stdio.h>
#include <iomanip>
#include "globals.h"

double adjConv=1e-7;
int adjMaxIt=960000;
int adjPrintIt=10000; 
int adjPrintConv=1; // 0 to hide real-time adjConvergence

void StegerJac(    std::vector <double> S,
                   std::vector <double> dx,
                   std::vector <double> W,
                   std::vector <double> &Ap_list,
                   std::vector <double> &An_list);

void adjFlux(      std::vector <double> S,
                   std::vector <double> dx,
                   std::vector <double> Ap_list,
                   std::vector <double> An_list,
                   std::vector <double> psi,
                   std::vector <double> &psiF);

void adjointEuler(  std::vector <double> S,
                    std::vector <double> dt,
                    std::vector <double> psiF,
                    std::vector <double> &Resi,
                    std::vector <double> &psi);

void adjointBC(std::vector <double> &psi, 
        std::vector <double> W,
        std::vector <double> dx,
        std::vector <double> S);
std::vector <double> adjoint(   std::vector <double> x, 
                                std::vector <double> dx, 
                                std::vector <double> S,
                                std::vector <double> W,
                                std::vector <double> &psi)
{
    std::vector <double> Resi(3*nx,0);
    std::vector <double> psiF(3*(nx+1),0);
    std::vector <double> dt(nx);

    std::vector <int> itV(adjMaxIt/adjPrintIt);
    std::vector <double> normV(adjMaxIt/adjPrintIt);

    std::vector <double> Ap_list(nx*3*3,0), An_list(nx*3*3,0);

    double normR=1.0;
    int iterations=0;


    adjointBC(psi, W, dx, S);
    StegerJac(S, dx, W, Ap_list, An_list);

    while(normR>adjConv && iterations<adjMaxIt)
    {
        iterations++;

        if(iterations%adjPrintIt==0) 
        {
            if(adjPrintConv==1)
            {
                std::cout<<"Iteration "<<iterations
                         <<"   NormR "<<std::setprecision(15)<<normR<<std::endl;
            }
            itV[iterations/adjPrintIt-1]=iterations;
            normV[iterations/adjPrintIt-1]=normR;
        }

        // CALCULATE TIME STEP
        for(int i=0;i<nx;i++)
            dt[i]=0.0001;


        adjFlux(S, dx, Ap_list, An_list, psi, psiF);
        adjointEuler(S, dt, psiF, Resi, psi);

        // Calculating the norm of the first costate residual
        normR=0;
        for(int i=0;i<nx;i++)
            normR=normR+Resi[2*nx+i]*Resi[2*nx+i];
        normR=sqrt(normR);
    }

    std::cout<<"Adjoint Iterations="<<iterations<<"   Density Residual="<<normR<<std::endl;
    std::vector <double> abc(1,0);
    
    FILE *Results;
    Results=fopen("Adjoint.dat","w");
    fprintf(Results,"%d\n",nx);
    for(int k=0;k<3;k++)
    for(int i=0;i<nx;i++)
        fprintf(Results, "%.15f\n",psi[k*nx+i]);

    fclose(Results);

    return abc;
}

void adjointEuler(std::vector <double> S,
                    std::vector <double> dt,
                    std::vector <double> psiF,
                    std::vector <double> &Resi,
                    std::vector <double> &psi)
{
    int ki;
    for(int k=0;k<3;k++)
    {
        for(int i=1;i<nx-1;i++)
        {
            ki=k*nx+i;
            Resi[ki] = psiF[ki+1] * S[i+1]
                       - psiF[ki] * S[i];
        }
        Resi[k*nx+0] = 0;
        Resi[k*nx+(nx-1)]= 0;
    }
    for(int k=0;k<3;k++)
    for(int i=1;i<nx-1;i++)
    {
        ki=k*nx+i;
        psi[ki] = psi[ki] - (dt[i]) * Resi[ki];
    }
    return;
}

void adjFlux(
                   std::vector <double> S,
                   std::vector <double> dx,
                   std::vector <double> Ap_list,
                   std::vector <double> An_list,
                   std::vector <double> psi,
                   std::vector <double> &psiF)
{
    // Calculate psiF
    for(int i=1; i<nx; i++)
    {
        psiF[0*nx+i]=0;
        psiF[1*nx+i]=0;
        psiF[2*nx+i]=0;
        for(int row=0;row<3;row++)
        {
            for(int col=0;col<3;col++)
            {
                int Ap_pos=((i-1)*3*3)+(row*3)+col;
                int An_pos=(i*3*3)+(row*3)+col;
                psiF[row*nx+i]=psiF[row*nx+i]
                                + Ap_list[Ap_pos]
                                * ( psi[col*nx+(i+1)] - psi[col*nx+i] )
                                + An_list[An_pos]
                                * ( psi[col*nx+i] - psi[col*nx+(i-1)] );
            }
        }
    }
}

// Calculates Jacobian
// Steger-Warming Flux Splitting
void StegerJac(    std::vector <double> S,
                   std::vector <double> dx,
                   std::vector <double> W,
                   std::vector <double> &Ap_list,
                   std::vector <double> &An_list)
{
    double eps=0.1;
    double gam=1.4;
    double M[3][3]={{0}},
           Minv[3][3]={{0}},
           N[3][3]={{0}},
           Ninv[3][3]={{0}},
           lambdaP[3][3],
           lambdaN[3][3];
    double lambdaa[3];
    
    
    double Ap[3][3], An[3][3], tempP[3][3], tempN[3][3], prefix[3][3], suffix[3][3];
    
    std::vector <double> rho(nx), u(nx), p(nx), c(nx);



    double beta=0.4;//gam-1;

    for(int i=0;i<nx;i++)
    {
        rho[i]=W[0*nx+i];
        u[i]=W[1*nx+i]/rho[i];
        p[i]=(gam-1)*(W[2*nx+i]-rho[i]*pow(u[i],2)/2);
        c[i]=sqrt((gam*p[i])/rho[i]);
    }


    for(int i=0;i<nx;i++)
    {
        for(int row=0;row<3;row++)
        for(int col=0;col<3;col++)
        {
            Ap[row][col]=0;
            An[row][col]=0;
            tempP[row][col]=0;
            tempN[row][col]=0;
            prefix[row][col]=0;
            suffix[row][col]=0;
            lambdaP[row][col]=0;
            lambdaN[row][col]=0;
        }
    
        M[0][0]=1;
        M[1][0]=-u[i]/rho[i];
        M[2][0]=0.5*u[i]*u[i]*beta;
        M[1][1]=1/rho[i];
        M[2][1]=-u[i]*beta;
        M[2][2]=beta;
        Minv[0][0]=1;
        Minv[1][0]=u[i];
        Minv[2][0]=0.5*u[i]*u[i];
        Minv[1][1]=rho[i];
        Minv[2][1]=u[i]*rho[i];
        Minv[2][2]=1/beta;
        N[0][0]=1;
        N[1][1]=rho[i]*c[i];
        N[2][1]=-rho[i]*c[i];
        N[0][2]=-1/(c[i]*c[i]);
        N[1][2]=1;
        N[2][2]=1;
        Ninv[0][0]=1;
        Ninv[0][1]=1/(2*c[i]*c[i]);
        Ninv[0][2]=1/(2*c[i]*c[i]);
        Ninv[1][1]=1/(2*rho[i]*c[i]);
        Ninv[1][2]=-1/(2*rho[i]*c[i]);
        Ninv[2][1]=0.5;
        Ninv[2][2]=0.5;
        lambdaa[0]=u[i];
        lambdaa[1]=u[i]+c[i];
        lambdaa[2]=u[i]-c[i];
        
        for(int k=0;k<3;k++)
            if(lambdaa[k]>0)
                lambdaP[k][k]=(lambdaa[k]+
                    sqrt(pow(lambdaa[k],2)+pow(eps,2)))/2;
            else
                lambdaN[k][k]=(lambdaa[k]-
                    sqrt(pow(lambdaa[k],2)+pow(eps,2)))/2;

        for(int row=0;row<3;row++)
        for(int col=0;col<3;col++)
            for(int k=0;k<3;k++)
            {
                prefix[row][col]+=Minv[row][k]*Ninv[k][col];
                suffix[row][col]+=N[row][k]*M[k][col];
            }
        for(int row=0;row<3;row++)
        for(int col=0;col<3;col++)
            for(int k=0;k<3;k++)
            {
                tempP[row][col]+=prefix[row][k]*lambdaP[k][col];
                tempN[row][col]+=prefix[row][k]*lambdaN[k][col];
            }
        for(int row=0;row<3;row++)
        for(int col=0;col<3;col++)
            for(int k=0;k<3;k++)
            {
                Ap[row][col]+=tempP[row][k]*suffix[k][col];
                An[row][col]+=tempN[row][k]*suffix[k][col];
            }
        // could remove above loop and just use aplist and anlist
        for(int row=0;row<3;row++)
        for(int col=0;col<3;col++)
        {
            int vec_pos=(i*3*3)+(row*3)+col;
            Ap_list[vec_pos]=Ap[row][col];
            An_list[vec_pos]=An[row][col];
        }

    }


}

// Set the adjoint BC for target pressure
void adjointBC(std::vector <double> &psi, 
        std::vector <double> W,
        std::vector <double> dx,
        std::vector <double> S)
{
    std::vector <double> pTarget(nx);

    ioTargetPressure(-1,pTarget);

    double gam=1.4;
    double rho0=W[0*nx];
    double u0=W[1*nx]/rho0;
    double p0=(gam-1)*(W[2*nx]-rho0*pow(u0,2)/2);
    double rhon=W[0*nx+nx-1];
    double un=W[1*nx+nx-1]/rhon;
    double pn=(gam-1)*(W[2*nx+nx-1]-rhon*pow(un,2)/2);

    psi[0*nx]=0;
    psi[1*nx-1]=0;
    psi[1*nx]= - (p0-pTarget[0])*dx[0]/(S[1]-S[0]);
    psi[2*nx-1]= - (pn-pTarget[nx-1])*dx[nx-1]/(S[nx]-S[nx-1]);
    psi[2*nx]=0;
    psi[3*nx-1]=0;
}
