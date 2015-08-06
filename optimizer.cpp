#include<iostream>
#include<math.h>
#include"quasiOneD.h"
#include"grid.h"
#include<vector>
#include<iomanip>

std::vector <double> finiteD(int nx, std::vector <double> x, std::vector <double> dx,
	       	std::vector <double> S, std::vector <double> designVar, int fitnessFun,
		double h, int method, double currentI);

void design(int nx, int descentType, int gradientType, int fitnessFun,
	      std::vector <double> x, std::vector <double> dx,
	       	std::vector <double> S, std::vector <double> designVar)
{
	int nDesVar=designVar.size();
	
	std::vector <double> normGradList;
	std::vector <double> B(nDesVar*nDesVar,0);
	double normGrad=1;
	double tolGrad=1e-8;
	int maxDesign=100;

	int printConv=1;

	double alpha=0.1;

	// 1 = Steepest Descent
	// 4 = QUASI-NEWTON (BFGS)
	
	// 1 = FD Forward
	// 2 = FD Backward
	// 3 = FD Centered

	double h=1e-3;
	std::vector <double> S_optimum;
	std::vector <double> gradient(nDesVar), pk(nDesVar), searchD(nDesVar);

	int rc;

	// Initialize B
	for(int r=0;r<nDesVar;r++)
	for(int c=0;c<nDesVar;c++)
	{
		if(r==c)
		{
			rc=r*nDesVar+c;
			B[rc]=1;
		}
		
	}

	normGradList.push_back(1);
	int iDesign=0;

	// Design Loop
	while(normGrad>tolGrad && iDesign<maxDesign)
	{
		iDesign++;


		if(printConv==1)
		{
			std::cout<<"Iteration :"<<iDesign<<
				"    GradientNorm: "<<normGrad<<std::endl;

			for(int i=0;i<nDesVar;i++)
				std::cout<<designVar[i]<<std::endl;


		}


		gradient=finiteD(nx,x,dx,S,designVar,fitnessFun,h,gradientType,-1);

		if(descentType==1)
		{
			for(int i=0;i<nDesVar;i++)
				pk[i]=-gradient[i];
		}

		for(int i=0;i<nDesVar;i++)
		{
			searchD[i]=alpha*pk[i];
		}

		for(int i=0;i<nDesVar;i++)
			designVar[i]=designVar[i]+searchD[i];

		normGrad=0;
		for(int i=0;i<nDesVar;i++)
			normGrad+=pow(gradient[i],2);
		normGrad=sqrt(normGrad);
		normGradList.push_back(normGrad);
		
	}

	for(int i=0;i<nDesVar;i++)
		std::cout<<gradient[i]<<std::endl;

	S=evalS(nx, designVar, x, dx);
	std::cout<<"Fitness: "<<quasiOneD(nx,x,dx,S,fitnessFun)<<std::endl;





	return;



}




std::vector <double> finiteD(int nx, std::vector <double> x, std::vector <double> dx,
	       	std::vector <double> S, std::vector <double> designVar, int fitnessFun,
	       	double h, int method, double currentI)
{
	// Method
	// 1 = Forward
	// 2 = Backward
	// 3 = Central
	std::vector <double> grad(nx);
	std::vector <double> tempS(nx+1);
	
	double I0, I1, I2, dh;

	int nDesVar=designVar.size();

	std::vector <double> tempD(nDesVar);

	if(currentI<0 && method!=3)
	{
		I0=quasiOneD(nx,x,dx,S,fitnessFun);
		std::cout<<"I0="<<std::setprecision(15)<<I0<<std::endl;
	}
	else
	{
		I0=currentI;	
	}

	for(int i=0;i<nDesVar;i++)
	{

		dh=designVar[i]*h;

		tempD=designVar;



		if(method==1)
		{
			tempD[i]+=dh;
			
			for(int k=0;k<nDesVar;k++)
				std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;

			tempS=evalS(nx, tempD, x, dx);

			I1=quasiOneD(nx,x,dx,tempS,fitnessFun);
			grad[i]=(I1-I0)/dh;
			std::cout<<"I1="<<std::setprecision(15)<<I1<<std::endl;
	
		}
		else if(method==2)
		{
			tempD[i]-=dh;
		
			for(int k=0;k<nDesVar;k++)
				std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;

			tempS=evalS(nx, tempD, x, dx);

			I2=quasiOneD(nx,x,dx,tempS,fitnessFun);
			grad[i]=(I0-I2)/dh;
			std::cout<<"I2="<<std::setprecision(15)<<I2<<std::endl;

		}
		else if(method==3)
		{
			tempD[i]+=dh;
			for(int k=0;k<nDesVar;k++)
				std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;


			tempS=evalS(nx, tempD, x, dx);


			I1=quasiOneD(nx,x,dx,tempS,fitnessFun);
			std::cout<<"I1="<<std::setprecision(15)<<I1<<std::endl;

			tempD=designVar;
			tempD[i]-=dh;

			for(int k=0;k<nDesVar;k++)
				std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;

			tempS=evalS(nx, tempD, x, dx);

			I2=quasiOneD(nx,x,dx,tempS,fitnessFun);
			std::cout<<"I2="<<std::setprecision(15)<<I2<<std::endl;

			grad[i]=(I1-I2)/(2*dh);
		}
		std::cout<<"Gradient "<<i+1<<":   "<<grad[i]<<std::endl<<std::endl;
	}

	return grad;
}

double stepBacktrackUncons(std::vector <double> designVar, std::vector <double> pk,
	       std::vector <double> gradient, double c1, double currentI,
	       int nx, std::vector <double> x, std::vector <double> dx,
	       int fitnessFun)
{
	double alpha=1;
	std::vector <double> tempS(nx+1);
	double newVal;

	double c_pk_grad=0;

	int nDesVar=designVar.size();

	for(int i=0;i<nDesVar;i++)
		c_pk_grad+=gradient[i]*pk[i];

	c_pk_grad*=c1;
	
	for(int i=0;nDesVar;i++)
		newVal=1;
	newVal=999999999999.0;


	for(int i=0;i<nDesVar;i++)
		tempD=designVar[i]+alpha*pk[i];
	tempS=evalS(nx,tempD,x,dx);
	newVal=quasiOneD(nx,x,dx,tempS,fitnessFun);

	while(newVal>(currentI+alpha*c_pk_grad))
	{
    		alpha=alpha*0.9;
	
		for(int i=0;i<nDesVar;i++)
			tempD=designVar[i]+alpha*pk[i];
		tempS=evalS(nx,tempD,x,dx);

	
	}

	return alpha;


	
}
