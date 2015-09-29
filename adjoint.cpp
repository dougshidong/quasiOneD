// Calculates the discrete costate fluxes
#include<iostream>
#include<math.h>
#include"quasiOneD.h"
#include<vector>
#include<iomanip>
#include<Eigen/Dense>


void adjointBC(int nx,
        std::vector <double> &psi, 
        std::vector <double> W,
        std::vector <double> dx,
        std::vector <double> S);

// Steger-Warming flux splitting scheme
std::vector <double> adjointSteger(int nx,
                                   std::vector <double> S,
                                   std::vector <double> dx,
                                   std::vector <double> W
                                   std::vector <double> psi,
                                   std::vector <double> &pFlux)

{
    double eps=0.1;
    double gam=1.4;
    double M[3][3]={0},
           Minv[3][3]={0},
           N[3][3]={0},
           Ninv[3][3]={0},
           lambdaP[3][3],
           lambdaN[3][3];
    double lambdaa[3];
    
    
    double Ap[3][3], An[3][3], tempP[3][3], tempN[3][3], prefix[3][3], suffix[3][3];
    
    std::vector <double> rho(nx), u(nx), p(nx), c(nx);

    std::vector <double> Ap_list(nx*3*3,0), An_list(nx*3*3,0);


    double beta=0.4;//gam-1;

    for(int i=1;i<nx-1;i++)
    {
        rho[i]=W[0*nx+i];       // rho
        u[i]=W[1*nx+i]/rho[i];  // U
        p[i]=(gam-1)*(W[2*nx+i]-rho[i]*pow(u[i],2)/2);  // Pressure
        c[i]=sqrt((gam*p[i])/rho[i]);// Speed of sound
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


    // Calculate Fluxes
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
                psiF[row*nx+i]=psiF[row*nx+i]+Ap_list[Ap_pos]*W[col*nx+(i-1)]
                               +An_list[An_pos]*W[col*nx+i];
            }
        }
    }
}


void adjointBC(int nx,
        std::vector <double> &psi, 
        std::vector <double> W,
        std::vector <double> dx,
        std::vector <double> S)
{
    std::vector <double> pTarget(nx);

    ioTargetPressure(-1,nx,pTarget);

    double gam=1.4;
    double rho0=W[0*nx];        // rho
    double u0=W[1*nx]/rho0; // U
    double p0=(gam-1)*(W[2*nx]-rho0*pow(u0,2)/2);  // Pressure
    double rhon=W[0*nx+nx-1];       // rho
    double un=W[1*nx+nx-1]/rhon;    // U
    double pn=(gam-1)*(W[2*nx+nx-1]-rhon*pow(un,2)/2);  // Pressure


    psi[2*nx+i]=(p0-pTarget[0])*dx[0]/(S[1]-S[0]);
    psi[nx-1]=(pn-pTarget[nx-1])*dx[nx-1]/(S[nx]-S[nx-1]);
}
