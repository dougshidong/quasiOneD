#include "quasiOneD.h"
#include "grid.h"
#include "optimizer.h"
#include <iostream>
#include <vector>

int main()
{
	//Geometry Parameters
	double a=0, b=1;
	double h=0.11, t1=0.5, t2=2.0;
//	double h=0.1, t1=0.7, t2=1.4;

	int nx=100;
	int fitnessFun=1;
	int descentType=4;

	int gradientType=1;

	std::vector <double> x(nx), S(nx+1);
	std::vector <double> dx(nx);
	std::vector <double> geom(3);
	std::vector <double> W(3*nx,0);
	

	geom[0]=h;
	geom[1]=t1;
	geom[2]=t2;
	
	x=evalX(nx,a,b);
	dx=evalDx(nx,x);
	S=evalS(nx,geom,x,dx);
/*
	std::cout<<"x\n";
	for(int i=0;i<x.size();i++)
		std::cout<<i<<" "<<x[i]<<std::endl;

	std::cout<<"dx\n";
	for(int i=0;i<dx.size();i++)
		std::cout<<i<<" "<<dx[i]<<std::endl;
	std::cout<<"S\n";
	for(int i=0;i<S.size();i++)
		std::cout<<i<<" "<<S[i]<<std::endl;
	std::cout<<"dx";

	for(int i=0;i<V.size();i++)
		std::cout<<i<<" "<<V[i]<<std::endl;
*/
	double fitness=quasiOneD(nx,x,dx,S,fitnessFun,geom, W);
	
//	design(nx, descentType, gradientType, fitnessFun, x, dx, S, geom);
	return 0;
}
