#include <iostream>
#include <math.h>
#include <Eigen/Dense>

// Discretization

const int nx=50;
double a=0, b=1;
double dx=(b-a)/(nx-1);

// Geometry

double h=0.15;
double t1=0.8;
double t2=3;

// Constants
double gam=1.4;
double R=286.9;
double Cv=R/(gam-1);

// Problem parameters

double Ttin=295.11;
double ptin=101192.6;
double pexit=0.8*ptin;

// Convergence Settings
double CFL=0.5;
double eps=0.3;
double normR=1;
double conv=1e-6;
int iterations=0;
int maxIt=20;



double isenP(double pt, double M);
double isenT(double Tt, double M);
void Flux_StegerWarming(double Flux[][nx-1], double W[][nx], double u[], double c[], double rho[]);

int quasiOneD()
{
	double x[nx],S[nx],V[nx];
	double rho[nx], u[nx], e[nx];
	double T[nx], p[nx], c[nx], Mach[nx];
	double W[3][nx], F[3][nx], Q[3][nx];

	double Flux[3][nx-1];
	double Resi[3][nx];
	double Resi1[3][nx],Resi2[3][nx],Resi3[3][nx];

	double dt[nx];
	double maxUC;


	// Initialize grid
	x[0]=a;
	x[nx-1]=b;
	for(int i=1; i<nx-1; i++)
		x[i]=dx/2+dx*(i-1);
		
	for(int i=0; i<nx; i++)
		S[i]=1-h*pow((sin(M_PI*pow(x[i],t1))),t2);

	for(int i=0; i<nx; i++)
		V[i]=(S[i]+S[i+1])/2 * dx;


	// Inlet flow properties
	Mach[0]=0.5;
	T[0]=isenT(Ttin,Mach[0]);
	p[0]=isenP(ptin,Mach[0]);
	rho[0]=p[0]/(R*T[0]);
	c[0]=sqrt(gam*p[0]/rho[0]);
	u[0]=Mach[0]*c[0];
	e[0]=rho[0]*(Cv*T[0]+0.5*pow(u[0],2));

	for(int i=1; i<nx; i++)
		Mach[i]=0.5;
	for(int i=1; i<nx; i++)
		p[i]=pexit;
	// Flow Properties Initialization
	for(int i=1; i<nx; i++)
	{
		T[i]=T[0];
		rho[i]=p[i]/(R*T[i]);
		c[i]=sqrt(gam*p[i]/rho[i]);
		u[i]=c[i]*Mach[i];
		e[i]=rho[i]*(Cv*T[i]+0.5*pow(u[i],2));
	}

	// State and Flux Vectors Initialization
	for(int i=0; i<nx; i++)
	{
		W[0][i]=rho[i];
		W[1][i]=rho[i]*u[i];
		W[2][i]=e[i];

		F[0][i]=rho[i]*u[i];
		F[1][i]=rho[i]*u[i]*u[i]+p[i];
		F[2][i]=(e[i]+p[i])*u[i];

		Q[0][i]=0;
		Q[2][i]=0;
	}
	Q[1][0]=0;
	Q[1][nx-1]=0;
	for(int i=1;i<nx-1;i++)
		Q[1][i]=p[i]*(S[i]-S[i-1]);


	while(normR>conv && iterations<maxIt)
	{
		iterations++;
		if(iterations%10==0) std::cout<<"Iteration "<<iterations<<std::endl;

		maxUC=0;
		for(int i=0;i<nx;i++)
			if(fabs(u[i]+c[i]>maxUC))
				maxUC=fabs(u[i]+c[i]);
		for(int i=0;i<nx;i++)
			dt[i]=(CFL*dx)/maxUC;

		Flux_StegerWarming(Flux,W,u,c,rho);

		// Euler Explicit
		//
		for(int k=0;k<3;k++)
		{
			for(int i=0;i<nx-1;i++)
				Resi[k][i+1]=Flux[k][i+1]*S[i+1]-Flux[k][i]*S[i]-Q[k][i+1];
			Resi[k][0]=0;
			Resi[k][nx-1]=0;
		}
		for(int k=0;k<3;k++)
			for(int i=0;i<nx-1;i++)
				W[k][i]=W[k][i]-(dt[i]/V[i-1])*Resi[k][i];
	}

	return 0;
}

double isenP(double pt, double M)
{
	return pt*pow((1+(gam-1)/2*pow(M,2)),(-gam/(gam-1)));
}

double isenT(double Tt, double M)
{
	return Tt*pow((1+(gam-1)/2*pow(M,2)),-1);
}


void matrixMult(double A[3][3], double B[3][3], double result[3][3])
{
    double temp[3][3];
    for(int row=0;row<3;row++)
        for(int col=0;col<3;col++)
            temp[row][col]=0;

    for(int row=0;row<3;row++)
        for(int col=0;col<3;col++)
        {
            for(int k=0;k<3;k++)
                temp[row][col]+=A[row][k]*B[k][col];
        }
    for(int row=0;row<3;row++)
        for(int col=0;col<3;col++)
            result[row][col]=temp[row][col];
}

void Flux_StegerWarming(double Flux[][nx-1], double W[][nx], double u[], double c[], double rho[])
{
	double S[nx-1][3][3],Sinv[nx-1][3][3],C[nx-1][3][3],Cinv[nx-1][3][3];
	double lambdaP[nx-1][3][3],lambdaN[nx-1][3][3];

	double Ap[nx-1][3][3],An[nx-1][3][3];
	double tempP[3][3], tempN[3][3], prefixM[3][3],suffixM[3][3];
	double beta=gam-1;
	double lambdaa[nx-1][3];


	memset(Flux,0,sizeof(Flux[0][0])*3*nx-1);
	memset(S,0,sizeof(S[0][0][0])*2*3*3);
	memset(Sinv,0,sizeof(Sinv[0][0][0])*2*3*3);
	memset(C,0,sizeof(C[0][0][0])*2*3*3);
	memset(Cinv,0,sizeof(Cinv[0][0][0])*2*3*3);
	memset(lambdaP,0,sizeof(lambdaP[0][0][0])*2*3*3);
	memset(lambdaN,0,sizeof(lambdaN[0][0][0])*2*3*3);


	for(int i=0;i<nx-1;i++)
	{
		S[i][0][0]=1;
		S[i][1][0]=-u[i]/rho[i];
		S[i][2][0]=0.5*u[i]*u[i]*beta;
		S[i][1][1]=1/rho[i];
		S[i][2][1]=-u[i]*beta;
		S[i][2][2]=beta;
		Sinv[i][0][0]=1;
		Sinv[i][1][0]=u[i];
		Sinv[i][2][0]=0.5*u[i]*u[i];
		Sinv[i][1][1]=rho[i];
		Sinv[i][2][1]=u[i]*rho[i];
		Sinv[i][2][2]=1/beta;
		C[i][0][0]=1;
		C[i][1][1]=rho[i]*c[i];
		C[i][2][1]=-rho[i]*c[i];
		C[i][0][2]=-1/(c[i]*c[i]);
		C[i][1][2]=1;
		C[i][2][2]=1;
		Cinv[i][0][0]=1;
		Cinv[i][0][1]=1/(2*c[i]*c[i]);
		Cinv[i][0][2]=1/(2*c[i]*c[i]);
		Cinv[i][1][1]=1/(2*rho[i]*c[i]);
		Cinv[i][1][2]=-1/(2*rho[i]*c[i]);
		Cinv[i][2][1]=0.5;
		Cinv[i][2][2]=0.5;
		lambdaa[i][0]=u[i];
		lambdaa[i][1]=u[i]+c[i];
		lambdaa[i][2]=u[i]-c[i];

		for(int k=0;k<3;k++)
			if(lambdaa[i][k]>0)
				lambdaP[i][k][k]=lambdaa[i][k]*
					sqrt(pow(lambdaa[i][k],2)+pow(eps,2))/2;
			else
				lambdaN[i][k][k]=lambdaa[i][k]*
					sqrt(pow(lambdaa[i][k],2)+pow(eps,2))/2;
	}

	for(int i=0;i<nx-1;i++)
	{
		memset(Ap,0,sizeof(Ap[0][0])*3*3);
		memset(An,0,sizeof(Ap[0][0])*3*3);
		for(int row=0;row<3;row++)
		for(int col=0;col<3;col++)
			for(int k=0;i<3;k++)
			{
				prefixM[row][col]+=Sinv[i][row][k]*Cinv[i][k][col];
				suffixM[row][col]+=C[i][row][k]*S[i][k][col];
			}
		for(int row=0;row<3;row++)
		for(int col=0;col<3;col++)
			for(int k=0;i<3;k++)
			{
				tempP[row][col]=prefixM[row][k]*lambdaP[i][k][col];
				tempN[row][col]=prefixM[row][k]*lambdaN[i][k][col];
			}
		for(int row=0;row<3;row++)
		for(int col=0;col<3;col++)
			for(int k=0;i<3;k++)
			{
				Ap[i][row][col]=tempP[row][k]*suffixM[k][col];
				An[i][row][col]=tempN[row][k]*suffixM[k][col];
			}
		for(int row=0;row<3;row++)
		for(int col=0;col<3;col++)
			Flux[row][i]=Flux[row][i]+Ap[i][row][col]*W[col][i]
				+An[i][row][col]*W[col][i+1];
	}


}


