#include <iostream>
#include <math.h>
#include <vector>

void InitializeGrid(int nx, std::vector <double> &x, std::vector <double> &dx, std::vector <double> &S)
{
	double h=0.15;
	double t1=0.8;
	double t2=3.0;
	double a=0, b=1;
	double dxConst=(b-a)/nx;
	
	for(int i=0; i<nx; i++)
		x.push_back(dxConst/2+dxConst*i);

	dx.push_back(x[1]-x[0]);

	for(int i=1; i<x.size()-1; i++)
	{
		dx.push_back( (x[i]-x[i-1])/2 + (x[i+1]-x[i])/2 );
	}
	dx.push_back(x[x.size()-1]-x[x.size()-2]);

	for(int i=0; i<x.size()-1; i++)
		S.push_back(1-h*pow(sin(M_PI*pow(x[i]-dx[i]/2,t1)),t2));
	
	S.push_back(1-h*pow(sin(M_PI*pow(x[x.size()]+dx[x.size()]/2,t1)),t2));

	return;
}
