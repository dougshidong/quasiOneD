#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>

const double PI=atan(1.0)*4.0;

// Constants
double gam=1.4;
double R=1716;
double Cv=R/(gam-1);

// Problem parameters

double Ttin=531.2;
double ptin=2117;
double pexit=0.72*ptin;
double a2=2*gam*Cv*Ttin*((gam-1)/(gam+1)); // used in isentropic nozzle


// Convergence Settings
double CFL=0.4;
double eps=0.3;
double normR=1.0;
double conv=1e-13;
int iterations=0;
int maxIt=200000;
int printIt=500;
int printConv=1;



double isenP(double pt, double M);
double isenT(double Tt, double M);
//void Flux_StegerWarming(double Flux[][nx-1], double W[][nx], double u[], double c[], double rho[]);

void Flux_StegerWarmingV(std::vector<std::vector<double> > &Flux, std::vector<std::vector<double> > W, std::vector<double> u, std::vector<double> c, std::vector<double> rho);



int quasiOneD(int nx, std::vector <double> x, std::vector <double> dx, std::vector <double> S,
		std::vector <double> V)
{
	std::vector <double> rho(nx), u(nx), e(nx);
	std::vector <double> T(nx), p(nx), c(nx), Mach(nx);
	std::vector <std::vector <double> > W(3,std::vector <double> (nx,0)),
					    F(3,std::vector <double> (nx,0)),
					    Q(3,std::vector <double> (nx,0));

	std::vector <std::vector <double> > Flux(3,std::vector <double> (nx+1));

	std::vector <std::vector <double> > Resi(3,std::vector <double> (nx,0)),
	  				    Resi1(3,std::vector <double> (nx,0)),
					    Resi2(3,std::vector <double> (nx,0)),
					    Resi3(3,std::vector <double> (nx,0));

	std::vector <double> dt(nx);
	double maxUC;

	double dpdu;
	std::vector <double> eigenvalues(3);
	std::vector <double> charRel(3);
	double dp, drho, du;
	double MachBound;


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
	if(iterations%printIt==0 && printConv==1) std::cout<<"Iteration "<<iterations<<"   NormR "<<std::setprecision(15)<<normR<<std::endl;

		maxUC=0;
		for(int i=0;i<nx;i++)
			if(fabs(u[i]+c[i]>maxUC))
				maxUC=fabs(u[i]+c[i]);
		for(int i=0;i<nx;i++)
			dt[i]=(CFL*dx[i])/maxUC;

		Flux_StegerWarmingV(Flux,W,u,c,rho);
//		for(int i=0;i<nx;i++)
//			std::cout<<Flux[1][i]<<std::endl;
/*
		for(int i=0;i<nx;i++)
		{
			std::cout<<"rho "<<rho[i]<<std::endl;
			std::cout<<"u "<<u[i]<<std::endl;
			std::cout<<"p "<<p[i]<<std::endl;
			std::cout<<"W "<<W[1][i]<<std::endl;
			std::cout<<"c "<<c[i]<<std::endl;
		}
*/


		// Euler Explicit
		//
		for(int k=0;k<3;k++)
		{
			for(int i=1;i<nx-1;i++)
				Resi[k][i]=Flux[k][i+1]*S[i+1]-Flux[k][i]*S[i]-Q[k][i];
			Resi[k][0]=0;
			Resi[k][nx-1]=0;
		}
		for(int k=0;k<3;k++)
			for(int i=0;i<nx-1;i++)
				W[k][i]=W[k][i]-(dt[i]/V[i-1])*Resi[k][i];
		// Inlet Boundary Condition
		//
		if(Mach[0]<1)
		{
			dpdu=ptin*(gam/(gam-1))
				*pow(1-((gam-1)/(gam+1))
				*u[0]*u[0]/a2,1/(gam-1))
				*(-2*((gam-1)/(gam+1))*u[0]/a2);
			eigenvalues[0]=((u[1]+u[0]-c[1]-c[0])/2)*(dt[0]/dx[0]);
			du=(-eigenvalues[0]*(p[1]-p[0]-rho[0]*c[0]*(u[1]-u[0])))
				/(dpdu-rho[0]*c[0]);

			u[0]=u[0]+du;
			T[0]=Ttin*(1-((gam-1)/(gam+1))*u[0]*u[0]/a2);
			p[0]=ptin*pow(T[0]/Ttin,gam/(gam-1));
			rho[0]=p[0]/(R*T[0]);
			e[0]=rho[0]*(Cv*T[0]+0.5*u[0]*u[0]);
			c[0]=sqrt(gam*p[0]/rho[0]);
			Mach[0]=u[0]/c[0];
		}

		// Exit boundary condition
		// NOTE NOT SURE IF DX IS FULL DX OR DX/2
		eigenvalues[0]=((u[nx-1]+u[nx-2])/2)*(dt[nx-1]/(dx[nx-1]));
		eigenvalues[1]=((u[nx-1]+u[nx-2])/2+(c[nx-1]+c[nx-2])/2)*(dt[nx-1]/(dx[nx-1]));
		eigenvalues[2]=((u[nx-1]+u[nx-2])/2-(c[nx-1]+c[nx-2])/2)*(dt[nx-1]/(dx[nx-1]));
		charRel[0]=-eigenvalues[0]*(rho[nx-1]-rho[nx-2]
				-(1/(c[nx-1]*c[nx-1]))*(p[nx-1]-p[nx-2]));
		charRel[1]=-eigenvalues[1]*(p[nx-1]-p[nx-2]+rho[nx-1]*c[nx-1]*(u[nx-1]-u[nx-2]));
		charRel[2]=-eigenvalues[2]*(p[nx-1]-p[nx-2]-rho[nx-1]*c[nx-1]*(u[nx-1]-u[nx-2]));

		MachBound=((u[nx-1]+u[nx-2])/2)/((c[nx-1]+c[nx-2])/2);
		//MachBound=(u[nx-1]/(c[nx-1]));
		if(MachBound>1)
		    dp=0.5*(charRel[1]+charRel[2]);
		else
		    dp=0;
	
		drho=charRel[0]+dp/(pow(c[nx-1],2));
		du=(charRel[1]-dp)/(rho[nx-1]*c[nx-1]);
	
		u[nx-1]=u[nx-1]+du;
		rho[nx-1]=rho[nx-1]+drho;
		p[nx-1]=p[nx-1]+dp;
		T[nx-1]=p[nx-1]/(rho[nx-1]*R);
		e[nx-1]=rho[nx-1]*(Cv*T[nx-1]+0.5*pow(u[nx-1],2));
		c[nx-1]=sqrt((gam*p[nx-1])/rho[nx-1]);
		Mach[nx-1]=u[nx-1]/c[nx-1];

		// Update flow properties
		//
		for(int i=1;i<nx-1;i++)
		{
		    rho[i]=W[0][i];				 // rho
		    u[i]=W[1][i]/W[0][i];			   // U
		    p[i]=(gam-1)*(W[2][i]-rho[i]*pow(u[i],2)/2);  // Pressure
		    T[i]=p[i]/(rho[i]*R);			   // Temperature
		    e[i]=W[2][i];				   // Energy
		    c[i]=sqrt((gam*p[i])/rho[i]);		 // Speed of sound
		    Mach[i]=u[i]/c[i];			      // Mach number
		}
		// Update vectors
		for(int i=0;i<nx;i++)
		{
		    W[0][i]=rho[i];
		    F[0][i]=rho[i]*u[i];
		    Q[0][i]=0;
		    W[1][i]=rho[i]*u[i];
		    F[1][i]=rho[i]*pow(u[i],2)+p[i];
		    W[2][i]=e[i];
		    F[2][i]=(e[i]+p[i])*u[i];
		    Q[2][i]=0;
		}
	
		for(int i=1;i<nx-1;i++)
		    Q[1][i]=p[i]*(S[i]-S[i-1]);
	
		// Calculating the norm of the density residual
		normR=0;
		for(int i=0;i<nx;i++)
		    normR=normR+Resi[0][i]*Resi[0][i];
		normR=sqrt(normR);
	}
	for(int k=0;k<3;k++)
	{
		std::cout<<"W"<<k+1<<std::endl;
		for(int i=0;i<nx;i++)
			std::cout<<W[k][i]<<std::endl;
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

/*
void Flux_StegerWarming(double Flux[][nx-1], double W[][nx], double u[], double c[], double rho[])
{
	double S[nx-1][3][3];
	double Sinv[nx-1][3][3];
	double C[nx-1][3][3];
	double Cinv[nx-1][3][3];
	double lambdaP[nx-1][3][3];
	double lambdaN[nx-1][3][3];

	double Ap[nx-1][3][3];
	double An[nx-1][3][3];
	double tempP[3][3];
	double tempN[3][3];
	double prefixMP[3][3];
	double suffixMP[3][3];
	double prefixMN[3][3];
	double suffixMN[3][3];

	double beta=gam-1;
	double lambdaa[nx-1][3];


	memset(Flux,0,sizeof(Flux[0][0])*3*(nx-1));
	memset(S,0,sizeof(S[0][0][0])*(nx-1)*3*3);
	memset(Sinv,0,sizeof(Sinv[0][0][0])*(nx-1)*3*3);
	memset(C,0,sizeof(C[0][0][0])*(nx-1)*3*3);
	memset(Cinv,0,sizeof(Cinv[0][0][0])*(nx-1)*3*3);
	memset(lambdaP,0,sizeof(lambdaP[0][0][0])*(nx-1)*3*3);
	memset(lambdaN,0,sizeof(lambdaN[0][0][0])*(nx-1)*3*3);
	memset(lambdaa,0,sizeof(lambdaa[0][0])*(nx-1)*3);

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
				lambdaP[i][k][k]=lambdaa[i][k]+
					sqrt(pow(lambdaa[i][k],2)+pow(eps,2))/2;
			else
				lambdaN[i][k][k]=lambdaa[i][k]-
					sqrt(pow(lambdaa[i][k],2)+pow(eps,2))/2;
	}

	for(int i=0;i<nx-1;i++)
	{
		memset(Ap,0,sizeof(Ap[0][0][0])*(nx-1)*3*3);
		memset(An,0,sizeof(An[0][0][0])*(nx-1)*3*3);
		memset(tempP,0,sizeof(tempP[0][0])*3*3);
		memset(tempN,0,sizeof(tempN[0][0])*3*3);
		memset(prefixMP,0,sizeof(prefixMP[0][0])*3*3);
		memset(suffixMP,0,sizeof(suffixMP[0][0])*3*3);
		memset(prefixMN,0,sizeof(prefixMN[0][0])*3*3);
		memset(suffixMN,0,sizeof(suffixMN[0][0])*3*3);

		for(int row=0;row<3;row++)
		for(int col=0;col<3;col++)
			for(int k=0;k<3;k++)
			{
				prefixMP[row][col]+=Sinv[i][row][k]*Cinv[i][k][col];
				suffixMP[row][col]+=C[i][row][k]*S[i][k][col];
				prefixMN[row][col]+=Sinv[i+1][row][k]*Cinv[i+1][k][col];
				suffixMN[row][col]+=C[i+1][row][k]*S[i+1][k][col];

			}
		for(int row=0;row<3;row++)
		for(int col=0;col<3;col++)
			for(int k=0;k<3;k++)
			{
				tempP[row][col]+=prefixMP[row][k]*lambdaP[i][k][col];
				tempN[row][col]+=prefixMN[row][k]*lambdaN[i+1][k][col];
			}
		for(int row=0;row<3;row++)
		for(int col=0;col<3;col++)
			for(int k=0;k<3;k++)
			{
				Ap[i][row][col]+=tempP[row][k]*suffixMP[k][col];
				An[i][row][col]+=tempN[row][k]*suffixMN[k][col];
			}
		for(int row=0;row<3;row++)
		for(int col=0;col<3;col++)
			Flux[row][i]=Flux[row][i]+Ap[i][row][col]*W[col][i]
				+An[i][row][col]*W[col][i+1];
	}


}
*/

void Flux_StegerWarmingV(std::vector<std::vector<double> > &Flux, std::vector<std::vector<double> > W, std::vector<double> u, std::vector<double> c, std::vector<double> rho)
{
	int nrho=rho.size();
	double S[3][3]={0},
	       Sinv[3][3]={0},
	       C[3][3]={0},
	       Cinv[3][3]={0},
	       lambdaP[3][3]={0},
	       lambdaN[3][3]={0};
	double lambdaa[3]={0};
	
	
	double Ap[3][3], An[3][3], tempP[3][3], tempN[3][3], prefix[3][3], suffix[3][3];
	


	std::vector <double> Ap_list(nrho*3*3,0), An_list(nrho*3*3,0);


	double beta=gam-1;


	for(int i=0;i<nrho;i++)
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
	
		S[0][0]=1;
		S[1][0]=-u[i]/rho[i];
		S[2][0]=0.5*u[i]*u[i]*beta;
		S[1][1]=1/rho[i];
		S[2][1]=-u[i]*beta;
		S[2][2]=beta;
		Sinv[0][0]=1;
		Sinv[1][0]=u[i];
		Sinv[2][0]=0.5*u[i]*u[i];
		Sinv[1][1]=rho[i];
		Sinv[2][1]=u[i]*rho[i];
		Sinv[2][2]=1/beta;
		C[0][0]=1;
		C[1][1]=rho[i]*c[i];
		C[2][1]=-rho[i]*c[i];
		C[0][2]=-1/(c[i]*c[i]);
		C[1][2]=1;
		C[2][2]=1;
		Cinv[0][0]=1;
		Cinv[0][1]=1/(2*c[i]*c[i]);
		Cinv[0][2]=1/(2*c[i]*c[i]);
		Cinv[1][1]=1/(2*rho[i]*c[i]);
		Cinv[1][2]=-1/(2*rho[i]*c[i]);
		Cinv[2][1]=0.5;
		Cinv[2][2]=0.5;
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
				prefix[row][col]+=Sinv[row][k]*Cinv[k][col];
				suffix[row][col]+=C[row][k]*S[k][col];
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

	for(int i=1; i<nrho; i++)
	{
		Flux[0][i]=0;
		Flux[1][i]=0;
		Flux[2][i]=0;
		for(int row=0;row<3;row++)
		for(int col=0;col<3;col++)
		{
			int Ap_pos=((i-1)*3*3)+(row*3)+col;
			int An_pos=(i*3*3)+(row*3)+col;
			Flux[row][i]=Flux[row][i]+Ap_list[Ap_pos]*W[col][i-1]
				+An_list[An_pos]*W[col][i];
		}
	}

}

