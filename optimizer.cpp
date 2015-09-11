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

std::vector <double> BFGS(int nDesVar, 
		std::vector <double> oldH, 
		std::vector <double> gradList, 
		std::vector <double> searchD);

void adjointBC(int nx,
		std::vector <double> &psi, 
		std::vector <double> W,
		std::vector <double> dx,
		std::vector <double> S);


void design(int nx, int descentType, int gradientType, int fitnessFun,
	      std::vector <double> x, std::vector <double> dx,
	       	std::vector <double> S, std::vector <double> designVar)
{
	std::vector <double> W(3*nx,0);

	int nDesVar=designVar.size();
	
	std::vector <double> normGradList, gradList;
	std::vector <double> H(nDesVar*nDesVar,0);
	double normGrad=1;
	double tolGrad=1e-5;
	double currentI;

	int maxDesign=10000;

	int printConv=1;

	double alpha=1;

	// 1 = Steepest Descent
	// 4 = QUASI-NEWTON (BFGS)
	
	// 1 = FD Forward
	// 2 = FD Backward
	// 3 = FD Centered

	double h=1e-4;
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
			std::cout<<"Current Design:\n";
			for(int i=0;i<nDesVar;i++)
				std::cout<<designVar[i]<<std::endl;


		}
		S=evalS(nx,designVar,x,dx);
//		currentI=quasiOneD(nx,x,dx,S,fitnessFun);
		currentI=quasiOneD(nx,x,dx,S,fitnessFun,designVar,W);

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
		std::cout<<"pk:\n";
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

	std::cout<<"Final Gradient:"<<std::endl;
	for(int i=0;i<nDesVar;i++)
		std::cout<<gradient[i]<<std::endl;

	std::cout<<std::endl<<"Final Design:"<<std::endl;
	for(int i=0;i<nDesVar;i++)
		std::cout<<designVar[i]<<std::endl;


	S=evalS(nx, designVar, x, dx);

	std::cout<<"Fitness: "<<quasiOneD(nx,x,dx,S,fitnessFun,designVar,W)<<std::endl;





	return;



}


std::vector <double> BFGS(int nDesVar, std::vector <double> oldH, std::vector <double> gradList, std::vector <double> searchD)
{
	int ls=gradList.size();
	int rc;
	std::vector <double> newH(nDesVar*nDesVar);
	double aa;
	Eigen::VectorXd dg(nDesVar), dx(nDesVar);
	Eigen::MatrixXd dH(nDesVar,nDesVar), a(nDesVar,nDesVar),b(nDesVar,nDesVar);

	for(int i=0;i<nDesVar;i++)
	{
		dg(i)=gradList[ls-nDesVar+i]-gradList[ls-2*nDesVar+i];
		dx(i)=searchD[i];
	}
	
	Eigen::Map <Eigen::Matrix<double,-1,-1,Eigen::RowMajor> > 
		cH(oldH.data(), nDesVar, nDesVar);

//	aa=( 1 + ((dg.transpose()*cH*dg) / (dg.transpose()*dx))(0) );
	
//	dH=aa*( (dx*dx.transpose())/(dx.transpose()*dg) )
//		- (cH*dg*dx.transpose() + (cH*dg*dx.transpose()).transpose())/(dg.transpose()*dx);


	a=((dx.transpose()*dg+dg.transpose()*cH*dg)(0)*(dx*dx.transpose()))
		/((dx.transpose()*dg)(0)*(dx.transpose()*dg)(0));
	b=(cH*dg*dx.transpose()+dx*dg.transpose()*cH)/(dx.transpose()*dg)(0);

	dH=a-b;

	std::cout<<"Current Inverse Hessian:"<<std::endl;
	for(int r=0;r<nDesVar;r++)
	{
		std::cout<<"\n";
		for(int c=0;c<nDesVar;c++)
		{
			rc=r*nDesVar+c;
			newH[rc]=oldH[rc]+dH(r,c);
			std::cout<<newH[rc]<<"\t\t";
		}
	}

	std::cout<<"\n\n";

	return newH;
}



void HessGradfiniteD(int nx, 
		std::vector <double> x, 
		std::vector <double> dx,
	       	std::vector <double> S, 
		std::vector <double> designVar, 
		int fitnessFun,
		double h, 
		int method, 
		double currentI, 
		std::vector <double> &G, 
		std::vector <double> H)
{
	std::vector <double> grad(nx);
	std::vector <double> tempS(nx+1);
	std::vector <double> W(3*nx,0);

	double I0, I1, dh1, dh2;

	int nDesVar=designVar.size();

	std::vector <double> tempD(nDesVar);

	if(currentI<0)
	{
//		I0=quasiOneD(nx,x,dx,S,fitnessFun);
		I0=quasiOneD(nx,x,dx,S,fitnessFun,designVar,W);

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
//			I1=quasiOneD(nx,x,dx,tempS,fitnessFun);
			I1=quasiOneD(nx,x,dx,tempS,fitnessFun,tempD,W);

			grad[i]=(I1-I0)/dh1;
			std::cout<<"I1="<<std::setprecision(15)<<I1<<std::endl;
		
		}
			std::cout<<"Gradient "<<i+1<<":   "<<grad[i]<<std::endl<<std::endl;
		}

}




std::vector <double> finiteD(int nx, std::vector <double> x, std::vector <double> dx,
	       	std::vector <double> S, std::vector <double> designVar, int fitnessFun,
	       	double h, int method, double currentI)
{
	std::vector <double> W(3*nx,0);

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
//		I0=quasiOneD(nx,x,dx,S,fitnessFun);
		I0=quasiOneD(nx,x,dx,S,fitnessFun,designVar,W);

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

//			I1=quasiOneD(nx,x,dx,tempS,fitnessFun);
			I1=quasiOneD(nx,x,dx,tempS,fitnessFun,tempD,W);

			grad[i]=(I1-I0)/dh;
			std::cout<<"I1="<<std::setprecision(15)<<I1<<std::endl;
	
		}
		else if(method==2)
		{
			tempD[i]-=dh;
		
			for(int k=0;k<nDesVar;k++)
				std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;

			tempS=evalS(nx, tempD, x, dx);

//			I2=quasiOneD(nx,x,dx,tempS,fitnessFun);
			I2=quasiOneD(nx,x,dx,tempS,fitnessFun,tempD,W);

			grad[i]=(I0-I2)/dh;
			std::cout<<"I2="<<std::setprecision(15)<<I2<<std::endl;

		}
		else if(method==3)
		{
			tempD[i]+=dh;
			for(int k=0;k<nDesVar;k++)
				std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;


			tempS=evalS(nx, tempD, x, dx);


//			I1=quasiOneD(nx,x,dx,tempS,fitnessFun);
			I1=quasiOneD(nx,x,dx,tempS,fitnessFun,tempD,W);

			std::cout<<"I1="<<std::setprecision(15)<<I1<<std::endl;

			tempD=designVar;
			tempD[i]-=dh;

			for(int k=0;k<nDesVar;k++)
				std::cout<<std::setprecision(15)<<tempD[k]<<std::endl;

			tempS=evalS(nx, tempD, x, dx);

//			I2=quasiOneD(nx,x,dx,tempS,fitnessFun);
			I2=quasiOneD(nx,x,dx,tempS,fitnessFun,tempD,W);

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
	std::vector <double> W(3*nx,0);

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
//	newVal=quasiOneD(nx,x,dx,tempS,fitnessFun);
	newVal=quasiOneD(nx,x,dx,tempS,fitnessFun,tempD,W);


	while(newVal>(currentI+alpha*c_pk_grad) && alpha>0.0001)
	{
    		alpha=alpha*0.5;
		std::cout<<"Alpha Reduction: "<<alpha<<std::endl;
	
		for(int i=0;i<nDesVar;i++)
			tempD[i]=designVar[i]+alpha*pk[i];
		tempS=evalS(nx,tempD,x,dx);
//		newVal=quasiOneD(nx,x,dx,tempS,fitnessFun);
		newVal=quasiOneD(nx,x,dx,tempS,fitnessFun,tempD,W);

	}

	return alpha;
}	


std::vector <double> adjointSteger(int nx,
				std::vector <double> S,
                                std::vector <double> dx,
				std::vector <double> W)
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
		rho[i]=W[0*nx+i];		// rho
		u[i]=W[1*nx+i]/rho[i];	// U
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

	std::vector <double> psi(nx);

        adjointBC(nx, psi, W, dx, S);

/*        
        for(int i=1; i<nx; i++)
        {
		Flux[0*nx+i]=0;
                Flux[1*nx+i]=0;
                Flux[2*nx+i]=0;
                for(int row=0;row<3;row++)
                for(int col=0;col<3;col++)
                {
                        int Ap_pos=((i-1)*3*3)+(row*3)+col;
                        int An_pos=(i*3*3)+(row*3)+col;
                        Flux[row*nx+i]=Flux[row*nx+i]+Ap_list[Ap_pos]*W[col*nx+(i-1)]
                                +An_list[An_pos]*W[col*nx+i];
                }
        }

*/


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
	double rho0=W[0*nx];		// rho
	double u0=W[1*nx]/rho0;	// U
	double p0=(gam-1)*(W[2*nx]-rho0*pow(u0,2)/2);  // Pressure
	double rhon=W[0*nx+nx-1];		// rho
	double un=W[1*nx+nx-1]/rhon;	// U
	double pn=(gam-1)*(W[2*nx+nx-1]-rhon*pow(un,2)/2);  // Pressure


	psi[0]=(p0-pTarget[0])*dx[0]/(S[1]-S[0]);
        psi[nx-1]=(pn-pTarget[nx-1])*dx[nx-1]/(S[nx]-S[nx-1]);
}


		
