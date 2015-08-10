#include<iostream>
#include<math.h>
#include"quasiOneD.h"
#include"grid.h"
#include<vector>
#include<iomanip>
#include<Eigen/Dense>

std::vector <double> finiteD(int nx, std::vector <double> x, std::vector <double> dx,
	       	std::vector <double> S, std::vector <double> designVar, int fitnessFun,
		double h, int method, double currentI);
double stepBacktrackUncons(std::vector <double> designVar, std::vector <double> pk,
	       std::vector <double> gradient, double currentI,
	       int nx, std::vector <double> x, std::vector <double> dx,
	       int fitnessFun);
std::vector <double> matrixMult(int mSize, std::vector <double> A, std::vector <double> B);

std::vector <double> BFGS(int nDesVar, 
		std::vector <double> oldH, 
		std::vector <double> gradList, 
		std::vector <double> searchD);


void design(int nx, int descentType, int gradientType, int fitnessFun,
	      std::vector <double> x, std::vector <double> dx,
	       	std::vector <double> S, std::vector <double> designVar)
{
	int nDesVar=designVar.size();
	
	std::vector <double> normGradList, gradList;
	std::vector <double> H(nDesVar*nDesVar,0);
	double normGrad=1;
	double tolGrad=1e-5;
	double currentI;

	int maxDesign=10;

	int printConv=1;

	double alpha=1;

	// 1 = Steepest Descent
	// 4 = QUASI-NEWTON (BFGS)
	
	// 1 = FD Forward
	// 2 = FD Backward
	// 3 = FD Centered

	double h=1e-7;
	std::vector <double> S_optimum;
	std::vector <double> gradient(nDesVar),
				pk(nDesVar),
				searchD(nDesVar),
				dgradient(nDesVar);

	int rc;
	
	// Initialize B
	for(int r=0;r<nDesVar;r++)
	for(int c=0;c<nDesVar;c++)
	{
		if(r==c)
		{
			rc=r*nDesVar+c;
			H[rc]=1;
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
		S=evalS(nx,designVar,x,dx);
		currentI=quasiOneD(nx,x,dx,S,fitnessFun);

		gradient=finiteD(nx,x,dx,S,designVar,fitnessFun,h,gradientType,currentI);
		
		for(int i=0;i<nDesVar;i++)
		{
			gradList.push_back(gradient[i]);
		}

		if(descentType==1)
		{
			for(int i=0;i<nDesVar;i++)
				pk[i]=-gradient[i];
		}
		if(descentType==4)
		{
			if(iDesign>1)
			{
				H=BFGS(nDesVar, H, gradList, searchD);
			}

			for(int r=0;r<nDesVar;r++)
			{
				pk[r]=0;
				for(int c=0;c<nDesVar;c++)
				{
					rc=r*nDesVar+c;
					pk[r]+=-H[rc]*gradient[c];
				}
			}
		}
		for(int i=0;i<nDesVar;i++)
			std::cout<<pk[i]<<std::endl;

                alpha=stepBacktrackUncons(designVar, pk, gradient, currentI,
					  nx, x, dx, fitnessFun);
		std::cout<<"Alpha="<<alpha<<std::endl;
		for(int i=0;i<nDesVar;i++)
		{
			searchD[i]=alpha*pk[i];
		}

		std::cout<<"Search Direction: :"<<std::endl;
		for(int i=0;i<nDesVar;i++)
			std::cout<<searchD[i]<<std::endl;




		for(int i=0;i<nDesVar;i++)
			designVar[i]=designVar[i]+searchD[i];

		normGrad=0;
		for(int i=0;i<nDesVar;i++)
			normGrad+=pow(gradient[i],2);
		normGrad=sqrt(normGrad);
		normGradList.push_back(normGrad);
		
		std::cout<<"End of Design Iteration: "<<iDesign<<std::endl<<std::endl<<std::endl;
	}

	for(int i=0;i<nDesVar;i++)
		std::cout<<gradient[i]<<std::endl;

	S=evalS(nx, designVar, x, dx);
	std::cout<<"Fitness: "<<quasiOneD(nx,x,dx,S,fitnessFun)<<std::endl;





	return;



}


std::vector <double> matrixMult(int mSize, std::vector <double> A, std::vector <double> B)
{
	double sum;
	int rc;
	int Ar=A.size()/mSize;
	int Bc=B.size()/mSize;

	std::vector <double> C(Ar*Bc,0);

	for(int r=0; r<Ar; r++)
	{
		for(int c=0; c<mSize; c++)
		{
			sum=0;
			rc=r*mSize+c;
			for(int k=0; k<mSize; k++)
			{
				C[rc]+=A[r*mSize+k]*B[k*mSize+c];
			}
		}
	}

	return C;
}

std::vector <double> BFGS(int nDesVar, std::vector <double> oldH, std::vector <double> gradList, std::vector <double> searchD)
{
	int ls=gradList.size();
	int rc;
	std::vector <double> newH(nDesVar*nDesVar);
	double aa;
	Eigen::VectorXd dg(nDesVar), dx(nDesVar);
	Eigen::MatrixXd dH(nDesVar,nDesVar);

	for(int i=0;i<nDesVar;i++)
	{
		dg(i)=gradList[(ls-1)-nDesVar+i]-gradList[(ls-1)-2*nDesVar+1];
		dx(i)=searchD[i];
	}
	
	//test=matrixMult(nDesVar,oldH,dg);

	Eigen::Map <Eigen::Matrix<double,-1,-1,Eigen::RowMajor> > 
		cH(oldH.data(), nDesVar, nDesVar);

	aa=( 1 + ((dg.transpose()*cH*dg) / (dg.transpose()*dx))(0) );
	
	dH=aa*( (dx*dx.transpose())/(dx.transpose()*dg) )
		- (cH*dg*dx.transpose() + (cH*dg*dx.transpose()).transpose())/(dg.transpose()*dx);

	for(int r=0;r<nDesVar;r++)
	for(int c=0;c<nDesVar;c++)
	{
		rc=r*nDesVar+c;
		newH[rc]=oldH[rc]+dH(r,c);
		std::cout<<newH[rc]<<std::endl;
	}

	return newH;
}



std::vector <double> HessianfiniteD(int nx, std::vector <double> x, std::vector <double> dx,
	       	std::vector <double> S, std::vector <double> designVar, int fitnessFun,
	       	double h, int method, double currentI)
{
	std::vector <double> grad(nx);
	std::vector <double> tempS(nx+1);
	
	double I0, I1, dh1, dh2;

	int nDesVar=designVar.size();

	std::vector <double> tempD(nDesVar);

	if(currentI<0)
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
		for(int j=i;j<nDesVar;j++)	
		{

			dh1=designVar[i]*h;
			dh2=designVar[j]*h;

			tempD=designVar;

			tempD[i]+=dh1;
			tempD[j]+=dh2;
			
			for(int k=0;k<nDesVar;k++)
				std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;
			
			tempS=evalS(nx, tempD, x, dx);
			I1=quasiOneD(nx,x,dx,tempS,fitnessFun);
			grad[i]=(I1-I0)/dh1;
			std::cout<<"I1="<<std::setprecision(15)<<I1<<std::endl;
		
		}
			std::cout<<"Gradient "<<i+1<<":   "<<grad[i]<<std::endl<<std::endl;
		}

	return grad;
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
	       std::vector <double> gradient, double currentI,
	       int nx, std::vector <double> x, std::vector <double> dx,
	       int fitnessFun)
{
	double alpha=1;
	double c1=1e-4;
	std::vector <double> tempS(nx+1);
	double newVal;

	double c_pk_grad=0;

	int nDesVar=designVar.size();
	std::vector <double> tempD(nDesVar);

	for(int i=0;i<nDesVar;i++)
	{
		c_pk_grad+=gradient[i]*pk[i];
	}

	c_pk_grad*=c1;
	
	for(int i=0;i<nDesVar;i++)
	{
		tempD[i]=designVar[i]+alpha*pk[i];
	}

	tempS=evalS(nx,tempD,x,dx);
	newVal=quasiOneD(nx,x,dx,tempS,fitnessFun);

	while(newVal>(currentI+alpha*c_pk_grad) && alpha>0.0001)
	{
    		alpha=alpha*0.5;
		std::cout<<"Alpha Reduction: "<<alpha<<std::endl;
	
		for(int i=0;i<nDesVar;i++)
			tempD[i]=designVar[i]+alpha*pk[i];
		tempS=evalS(nx,tempD,x,dx);
		newVal=quasiOneD(nx,x,dx,tempS,fitnessFun);
	}

	return alpha;


	
}
