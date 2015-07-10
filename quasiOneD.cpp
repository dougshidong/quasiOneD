#include <iostream>
#include <math.h>
#include <Eigen/Dense>

// Discretization

int nx=50;
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
double normR=1;
double conv=1e-6;
double iterations=0;
double maxIt=20;



double isenP(double pt, double M);
double isenT(double Tt, double M);

int quasiOneD()
{
	double x[nx],S[nx],V[nx];
	double rho[nx], u[nx], e[nx];
	double T[nx], p[nx], c[nx], Mach[nx];
	double W[3][nx], F[3][nx], Q[3][nx];

	double dt[nx];
	double maxUc;


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

		maxUc=0;
		for(int i=0;i<nx;i++)
			if(fabs(u[i]+c[i]>maxUc))
				maxUC=fabs(u[i]+c[i];
		for(int i=0;i<nx;i++)
			dt[i]=(CFL*dx)/maxUc;

		Flux_StegerWarming(Flux,F,U,C,W);
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

//void scalarF(double **Flux, double **F, double U[], double c[], double W[][])
//{
//}


