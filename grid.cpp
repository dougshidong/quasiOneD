#include <iostream>
#include <math.h>
#include <vector>
#include "globals.h"

// Evaluate X
std::vector <double> evalX(double a, double b)
{
	std::vector <double> x(nx);
	double dxConst=(b-a)/nx;
	
	for(int i=0; i<nx; i++)
		x[i]=dxConst/2+dxConst*i;

	return x;
}

//  Evaluate dx

std::vector <double> evalDx(std::vector <double> x)
{
	std::vector <double> dx(nx);
	
	dx[0]=x[1]-x[0];
	for(int i=1; i<nx; i++)
	{
		dx[i]= (x[i]-x[i-1])/2 + (x[i+1]-x[i])/2 ;
	}
	dx[nx-1]=x[x.size()-1]-x[x.size()-2];

	return dx;
}

std::vector <double> evalS(std::vector <double> geom,
			std::vector <double> x, std::vector <double> dx)
{
	std::vector <double> S(nx+1);

	// Define Area
	for(int i=0; i<nx; i++)
		S[i]= 1-geom[0]*pow(sin(M_PI*pow(fabs(x[i]-dx[i]/2),geom[1])),geom[2]);
	
	S[nx]= 1-geom[0]*pow(sin(M_PI*pow(x[nx-1]+dx[nx-1]/2,geom[1])),geom[2]);

	return S;


}

void InitializeGrid(std::vector <double> &x, std::vector <double> &dx,
		std::vector <double> &S)
{
	double h=0.15;
	double t1=0.8;
	double t2=3.0;
	double a=0, b=1;
	double dxConst=(b-a)/nx;
	
	// Define X
	for(int i=0; i<nx; i++)
		x.push_back(dxConst/2+dxConst*i);

	// Define dx
	dx.push_back(x[1]-x[0]);
	for(int i=1; i<nx-1; i++)
	{
		dx.push_back( (x[i]-x[i-1])/2 + (x[i+1]-x[i])/2 );
	}
	dx.push_back(x[x.size()-1]-x[x.size()-2]);

	// Define Area
	for(int i=0; i<nx; i++)
		S.push_back(1-h*pow(sin(M_PI*pow(fabs(x[i]-dx[i]/2),t1)),t2));
	
	S.push_back(1-h*pow(sin(M_PI*pow(x[x.size()]+dx[x.size()]/2,t1)),t2));

	return;
}
