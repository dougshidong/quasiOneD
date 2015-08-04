#include<stdio.h>
#include<math.h>
#include<iostream>
const double PI = atan(1.0)*4;

// Problem definition
double h=0.15;
double t1=0.8;
double t2=3.0;

double gam=1.4;
double Tt=531.2;
double pt=2117;
double R=1716;
double Minlet=1.2;
double pexit=0.72*pt;
double Cv=R/(gam-1);
double a2=2*gam*Cv*Tt*((gam-1)/(gam+1)); // used in isentropic nozzle

double a=0,b=1;         // Bounds on x
const int nx=50;        // Number of grid points

double eps=0.3; // Epsilon
double CFL=0.4; // CFL
int maxIt=200000;

// Geometry areas
double shape(double x)
{
    return 1-h*(pow(sin(PI*pow(x,t1)),t2));
}

// Subsonic inlet boundary
double SubIsenInP(double m)
{
    return pt*pow(1+m*m*(gam-1)/2,-gam/(gam-1));
}
double SubIsenInT(double m)
{
    return Tt/(1+m*m*(gam-1)/2);
}

double max3NonAbs(double A, double B, double C)
{
    if((A)>(B))
    {
        if((A)>(C))
            return A;
        else
            return C;
    }
    else
    {
        if((B)>(C))
            return B;
    }
    return C;
}
double max3(double A, double B, double C)
{
    if(fabs(A)>fabs(B))
    {
        if(fabs(A)>fabs(C))
            return fabs(A);
        else
            return fabs(C);
    }
    else
    {
        if(fabs(B)>fabs(C))
            return fabs(B);
    }
    return fabs(C);
}

void scalarF(double Flux[][nx-1], double F[][nx], double U[], double c[], double W[][nx])
{
    for(int k=0;k<3;k++)
        for(int i=0;i<nx-1;i++)
        {
            Flux[k][i]=0.5*(F[k][i]+F[k][i+1])-0.5*eps*
                max3((U[i]+U[i+1])/2,((U[i]+U[i+1])/2+(c[i]+c[i+1])/2),((U[i]+U[i+1])/2-(c[i]+c[i+1])/2))*
                (W[k][i+1]-W[k][i]);
        }
}

void matrixMult(double A[][3], double B[][3], double result[][3])
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

void StegerWarmingF(double Flux[][nx-1], double W[][nx],double U[], double c[],double rho[])
{
    double S[3][3], Sinv[3][3], C[3][3],Cinv[3][3], lambdaP[3][3],lambdaN[3][3];
    double S2[3][3], Sinv2[3][3], C2[3][3],Cinv2[3][3], lambdaP2[3][3],lambdaN2[3][3];
    double Ap[3][3],An[3][3];
    double beta=gam-1;
    double eigenvalue[3];
    double eigenvalue2[3];
    for(int i=0;i<nx-1;i++)
        for(int k=0;k<3;k++)
            Flux[k][i]=0;

    for(int i=0;i<nx-1;i++)
    {
        for(int row=0;row<3;row++)
            for(int col=0;col<3;col++)
            {
                S[row][col]=0;
                Sinv[row][col]=0;
                C[row][col]=0;
                Cinv[row][col]=0;
                lambdaP[row][col]=0;
                lambdaN[row][col]=0;
                Ap[row][col]=0;
                An[row][col]=0;
                S2[row][col]=0;
                Sinv2[row][col]=0;
                C2[row][col]=0;
                Cinv2[row][col]=0;
                lambdaP2[row][col]=0;
                lambdaN2[row][col]=0;
            }
        S[0][0]=1;
        S[1][0]=-U[i]/rho[i];
        S[2][0]=0.5*U[i]*U[i]*beta;
        S[1][1]=1/rho[i];
        S[2][1]=-U[i]*beta;
        S[2][2]=beta;
        Sinv[0][0]=1;
        Sinv[1][0]=U[i];
        Sinv[2][0]=0.5*U[i]*U[i];
        Sinv[1][1]=rho[i];
        Sinv[2][1]=U[i]*rho[i];
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
        eigenvalue[0]=U[i];
        eigenvalue[1]=U[i]+c[i];
        eigenvalue[2]=U[i]-c[i];
        for(int k=0;k<3;k++)
        {
            if(eigenvalue[k]>0)
                lambdaP[k][k]=(eigenvalue[k]+sqrt(pow(eigenvalue[k],2)+pow(eps,2)))/2;
            else
                lambdaN[k][k]=(eigenvalue[k]-sqrt(pow(eigenvalue[k],2)+pow(eps,2)))/2;
        }
        S2[0][0]=1;
        S2[1][0]=-U[i+1]/rho[i+1];
        S2[2][0]=0.5*U[i+1]*U[i+1]*beta;
        S2[1][1]=1/rho[i+1];
        S2[2][1]=-U[i+1]*beta;
        S2[2][2]=beta;
        Sinv2[0][0]=1;
        Sinv2[1][0]=U[i+1];
        Sinv2[2][0]=0.5*U[i+1]*U[i+1];
        Sinv2[1][1]=rho[i+1];
        Sinv2[2][1]=U[i+1]*rho[i+1];
        Sinv2[2][2]=1/beta;
        C2[0][0]=1;
        C2[1][1]=rho[i+1]*c[i+1];
        C2[2][1]=-rho[i+1]*c[i+1];
        C2[0][2]=-1/(c[i+1]*c[i+1]);
        C2[1][2]=1;
        C2[2][2]=1;
        Cinv2[0][0]=1;
        Cinv2[0][1]=1/(2*c[i+1]*c[i+1]);
        Cinv2[0][2]=1/(2*c[i+1]*c[i+1]);
        Cinv2[1][1]=1/(2*rho[i+1]*c[i+1]);
        Cinv2[1][2]=-1/(2*rho[i+1]*c[i+1]);
        Cinv2[2][1]=0.5;
        Cinv2[2][2]=0.5;
        eigenvalue2[0]=U[i+1];
        eigenvalue2[1]=U[i+1]+c[i+1];
        eigenvalue2[2]=U[i+1]-c[i+1];
        for(int k=0;k<3;k++)
        {
            if(eigenvalue2[k]>0)
                lambdaP2[k][k]=(eigenvalue2[k]+sqrt(pow(eigenvalue2[k],2)+pow(eps,2)))/2;
            else
                lambdaN2[k][k]=(eigenvalue2[k]-sqrt(pow(eigenvalue2[k],2)+pow(eps,2)))/2;
        }

        matrixMult(Sinv,Cinv,Ap);
        matrixMult(Ap,lambdaP,Ap);
        matrixMult(Ap,C,Ap);
        matrixMult(Ap,S,Ap);

        matrixMult(Sinv2,Cinv2,An);
        matrixMult(An,lambdaN2,An);
        matrixMult(An,C2,An);
        matrixMult(An,S2,An);

        for(int row=0;row<3;row++)
        {
            for(int col=0;col<3;col++)
            {
                Flux[row][i]=Flux[row][i]+Ap[row][col]*W[col][i]+An[row][col]*W[col][i+1];
            }
        }
    }

}

void modifiedStegerWarmingF(double Flux[][nx-1], double W[][nx],double U[], double c[],double rho[])
{
    double S[3][3], Sinv[3][3], C[3][3],Cinv[3][3], lambdaP[3][3],lambdaN[3][3];
    double Ap[3][3],An[3][3];
    double beta=gam-1;
    double eigenvalue[3];
    for(int i=0;i<nx-1;i++)
        for(int k=0;k<3;k++)
            Flux[k][i]=0;

    for(int i=0;i<nx-1;i++)
    {
        for(int row=0;row<3;row++)
            for(int col=0;col<3;col++)
            {
                S[row][col]=0;
                Sinv[row][col]=0;
                C[row][col]=0;
                Cinv[row][col]=0;
                lambdaP[row][col]=0;
                lambdaN[row][col]=0;
                Ap[row][col]=0;
                An[row][col]=0;
            }
        S[0][0]=1;
        S[1][0]=(-U[i]/rho[i]-U[i+1]/rho[i+1])/2;
        S[2][0]=(0.5*U[i]*U[i]*beta+0.5*U[i+1]*U[i+1]*beta)/2;
        S[1][1]=(1/rho[i]+1/rho[i+1])/2;
        S[2][1]=(-U[i]*beta-U[i+1]*beta)/2;
        S[2][2]=beta;
        Sinv[0][0]=1;
        Sinv[1][0]=(U[i]+U[i+1])/2;
        Sinv[2][0]=(0.5*U[i]*U[i]+0.5*U[i+1]*U[i+1])/2;
        Sinv[1][1]=(rho[i]+rho[i+1])/2;
        Sinv[2][1]=(U[i]*rho[i]+U[i+1]*rho[i+1])/2;
        Sinv[2][2]=1/beta;
        C[0][0]=1;
        C[1][1]=(rho[i]*c[i]+rho[i+1]*c[i+1])/2;
        C[2][1]=(-rho[i]*c[i]-rho[i+1]*c[i+1])/2;
        C[0][2]=(-1/(c[i]*c[i])-1/(c[i+1]*c[i+1]))/2;
        C[1][2]=1;
        C[2][2]=1;
        Cinv[0][0]=1;
        Cinv[0][1]=(1/(2*c[i]*c[i])+1/(2*c[i+1]*c[i+1]))/2;
        Cinv[0][2]=(1/(2*c[i]*c[i])+1/(2*c[i+1]*c[i+1]))/2;
        Cinv[1][1]=(1/(2*rho[i]*c[i])+1/(2*rho[i+1]*c[i+1]))/2;
        Cinv[1][2]=(-1/(2*rho[i]*c[i])-1/(2*rho[i+1]*c[i+1]))/2;
        Cinv[2][1]=0.5;
        Cinv[2][2]=0.5;
        eigenvalue[0]=(U[i]+U[i+1])/2;
        eigenvalue[1]=(U[i]+c[i]+U[i+1]+c[i+1])/2;
        eigenvalue[2]=(U[i]-c[i]+U[i+1]-c[i+1])/2;
        for(int k=0;k<3;k++)
        {
            if(eigenvalue[k]>0)
                lambdaP[k][k]=(eigenvalue[k]+sqrt(pow(eigenvalue[k],2)+pow(eps,2)))/2;
            else
                lambdaN[k][k]=(eigenvalue[k]-sqrt(pow(eigenvalue[k],2)+pow(eps,2)))/2;
        }

        matrixMult(Sinv,Cinv,Ap);
        matrixMult(Ap,lambdaP,Ap);
        matrixMult(Ap,C,Ap);
        matrixMult(Ap,S,Ap);

        matrixMult(Sinv,Cinv,An);
        matrixMult(An,lambdaN,An);
        matrixMult(An,C,An);
        matrixMult(An,S,An);

        for(int row=0;row<3;row++)
        {
            for(int col=0;col<3;col++)
            {
                Flux[row][i]=Flux[row][i]+Ap[row][col]*W[col][i]+An[row][col]*W[col][i+1];
            }
        }
    }
}

double minn(double a, double b)
{
    if(fabs(a)>fabs(b))
        return a;
    return b;
}
void correctedModifiedStegerWarmingF(double Flux[][nx-1], double W[][nx],double U[], double c[],double rho[],double p[])
{
    double MSW[3][nx-1], SW[3][nx-1];
    double Whalf[3][nx];
    modifiedStegerWarmingF(MSW,W,U,c,rho);
    StegerWarmingF(SW,W,U,c,rho);
    for(int k=0;k<3;k++)
        for(int i=0;i<nx-1;i++)
        {
            Whalf[k][i]=1/(1+pow((p[i+1]-p[i])/minn(p[i+1],p[i]),2));
        }
    for(int k=0;k<3;k++)
        for(int i=0;i<nx-1;i++)
            Flux[k][i]=Whalf[k][i]*MSW[k][i]+(1-Whalf[k][i])*SW[k][i];
}

void RoeF(double Flux[][nx-1], double W[][nx],double F[][nx], double U[], double c[],double rho[],double p[])
{
    double epsilon;
    double temp;
    double rhoH;
    double uH;
    double hH;
    double cH;
    double A[3][3];
    double S[3][3], Sinv[3][3], C[3][3],Cinv[3][3], lambda[3][3];
    double beta=gam-1;
    double eigenvalueH[3];
    double eigenvalueI[3];
    double eigenvalueIp1[3];

    for(int i=0;i<nx-1;i++)
        for(int k=0;k<3;k++)
            Flux[k][i]=0;

    for(int i=0;i<nx-1;i++)
    {
        rhoH=sqrt(rho[i]*rho[i+1]);
        uH=(sqrt(rho[i])*U[i]+sqrt(rho[i+1])*U[i+1])/(sqrt(rho[i])+sqrt(rho[i+1]));
        hH=(sqrt(rho[i])*((W[2][i]+p[i])/rho[i])+sqrt(rho[i+1])*((W[2][i+1]+p[i+1])/rho[i+1]))/(sqrt(rho[i])+sqrt(rho[i+1]));
        cH=sqrt((gam-1)*(hH-0.5*uH*uH));
    }

    for(int i=0;i<nx-1;i++)
    {
        for(int row=0;row<3;row++)
            for(int col=0;col<3;col++)
            {
                S[row][col]=0;
                Sinv[row][col]=0;
                C[row][col]=0;
                Cinv[row][col]=0;
                lambda[row][col]=0;
                A[row][col]=0;
            }
        S[0][0]=1;
        S[1][0]=-uH/rhoH;
        S[2][0]=0.5*uH*uH*beta;
        S[1][1]=1/rhoH;
        S[2][1]=-uH*beta;
        S[2][2]=beta;
        Sinv[0][0]=1;
        Sinv[1][0]=uH;
        Sinv[2][0]=0.5*uH*uH;
        Sinv[1][1]=rhoH;
        Sinv[2][1]=uH*rhoH;
        Sinv[2][2]=1/beta;
        /*
        C[0][0]=1;
        C[1][1]=rhoH*cH;
        C[2][1]=-rhoH*cH;
        C[0][2]=-1/(cH*cH);
        C[1][2]=1;
        C[2][2]=1;
        Cinv[0][0]=1;
        Cinv[0][1]=1/(2*cH*cH);
        Cinv[0][2]=1/(2*cH*cH);
        Cinv[1][1]=1/(2*rhoH*cH);
        Cinv[1][2]=-1/(2*rhoH*cH);
        Cinv[2][1]=0.5;
        Cinv[2][2]=0.5;*/
        eigenvalueH[0]=uH;
        eigenvalueH[1]=uH+cH;
        eigenvalueH[2]=uH-cH;
        eigenvalueI[0]=U[i];
        eigenvalueI[1]=U[i]+c[i];
        eigenvalueI[2]=U[i]-c[i];
        eigenvalueIp1[0]=U[i+1];
        eigenvalueIp1[1]=U[i+1]+c[i+1];
        eigenvalueIp1[2]=U[i+1]-c[i+1];

        for(int k=0;k<3;k++)
        {
            epsilon=max3NonAbs(0,eigenvalueH[k]-eigenvalueI[k],eigenvalueIp1[k]-eigenvalueH[k]);
            if(fabs(eigenvalueH[k])<=epsilon)
                eigenvalueH[k]=0.5*(eigenvalueH[k]/epsilon+epsilon);
        }

        for(int k=0;k<3;k++)
        {
            lambda[k][k]=eigenvalueH[k];
        }


        /*
        matrixMult(Sinv,Cinv,A);
        matrixMult(A,lambda,A);
        matrixMult(A,C,A);
        matrixMult(A,S,A);
        */
        matrixMult(Sinv,lambda,A);
        matrixMult(A,S,A);

        for(int k=0;k<3;k++)
            for(int i=0;i<nx-1;i++)
            {
                temp=0;
                for(int j=0;j<3;j++)
                {
                    temp+=fabs(A[k][j])*(W[j][i+1]-W[j][i]);
                }
                Flux[k][i]=0.5*(F[k][i]+F[k][i+1])-0.5*temp;
            }
    }



}

int main()
{
    FILE *xData, *DensityResiNorm, *PressureDist, *MachDist, *IterationsData, *PressureLoss, *EnergyDist;
    DensityResiNorm=fopen("DensityResiNorm.csv","w");
    PressureDist=fopen("PressureDist.csv","w");
    xData=fopen("xData.csv","w");
    MachDist=fopen("MachDist.csv","w");
    EnergyDist=fopen("EnergyDist.csv","w");
    IterationsData=fopen("IterationsData.csv","w");
    PressureLoss=fopen("PressureLoss.csv","w");
    FILE *test, *test2;
    test=fopen("test.csv","w");
    test2=fopen("test2.csv","w");
    double dx=(b-a)/(nx-2);
    double x[nx];
    double S[nx-1];
    double V[nx-2];
    double T[nx], e[nx], U[nx], p[nx], rho[nx], c[nx], Mach[nx];
    double W[3][nx], F[3][nx], Q[3][nx];

    double W0[3][nx],W1[3][nx],W2[3][nx],W3[3][nx];
    double F0[3][nx],F1[3][nx],F2[3][nx];
    double Q0[3][nx],Q1[3][nx],Q2[3][nx];

    double Utemp[nx];
    double ptemp[nx], rhotemp[nx], ctemp[nx];

    double Flux[3][nx-1];
    double Resi[3][nx];
    double Resi0[3][nx],Resi1[3][nx],Resi2[3][nx];

    double dpdu;
    double eigenvalues[3];
    double charRel[3];
    double MachBound;
    double dp, drho, du;
    double timestep;
    double maxUc;
    double dt[nx];
    double normm=1;
    double conv=pow(10,-13);
    int iterations=0;

    // Discretization
    x[0]=a;
    x[nx-1]=b;
    for(int i=1;i<nx-1;i++)
        x[i]=(dx/2)+(i-1)*dx;
    for(int i=0;i<nx;i++)
        fprintf(xData,"%f\n",x[i]);
    // Area at surface of each CV
    for(int i=0;i<nx-1;i++)
        S[i]=shape(i*dx);
    // Control volumes
    for(int i=0;i<nx-2;i++)
        V[i]=dx*(S[i]+S[i+1])/2;

    // Inlet flow properties
    Mach[0]=Minlet;
    T[0]=SubIsenInT(Mach[0]);
    p[0]=SubIsenInP(Mach[0]);
    rho[0]=p[0]/(R*T[0]);
    c[0]=sqrt((gam*p[0])/rho[0]);
    U[0]=Mach[0]*c[0];
    e[0]=rho[0]*(Cv*T[0]+0.5*pow(U[0],2));

    // Velocity Initialization
    for(int i=1;i<nx;i++)
        Mach[i]=Minlet;
    for(int i=1;i<nx;i++)
        p[i]=pexit;//SubIsenInP(Mach[i]);                              // Pressure
    p[nx-1]=pexit;
    //Mach[nx-1]=0.9;
    // Flow properties inititialization     // Assuming exit has constant pressure
    for(int i=1;i<nx;i++)
    {
        //T[i]=SubIsenInT(Mach[i]);
        T[i]=T[0];                              // Temperature
        rho[i]=p[i]/(R*T[i]);                   // rho
        c[i]=sqrt((gam*p[i])/rho[i]);         // Speed of sound
        U[i]=c[i]*Mach[i];
        e[i]=rho[i]*(Cv*T[i]+0.5*pow(U[i],2));  // Energy
    }

    for(int i=0;i<nx;i++)
    {
        W[0][i]=rho[i];
        F[0][i]=rho[i]*U[i];
        Q[0][i]=0;
        W[1][i]=rho[i]*U[i];
        F[1][i]=rho[i]*pow(U[i],2)+p[i];
        W[2][i]=e[i];
        F[2][i]=(e[i]+p[i])*U[i];
        Q[2][i]=0;
    }
    Q[1][0]=0;
    Q[1][nx-1]=0;
    for(int i=1;i<nx-1;i++)
        Q[1][i]=p[i]*(S[i]-S[i-1]);


    for(int i=0;i<nx;i++)
        fprintf(test,"%f,",Mach[i]);
    fprintf(test,"\n");

    for(int i=0;i<nx;i++)
        fprintf(test2,"%f,",p[i]/pt);
    fprintf(test2,"\n");

    // ITERATIONS
    while(normm>conv && iterations<maxIt)
    {
        iterations++;

        if(iterations%500==0)	std::cout<<"iteration"<<iterations<<" normm: "<<normm<<std::endl;
        fprintf(IterationsData,"%d\n",iterations);
        maxUc=0;
        for(int i=0;i<nx;i++)
            if(fabs(U[i]+c[i]>maxUc))
               maxUc=fabs(U[i]+c[i]);
        // Time step for each CV
        for(int i=0;i<nx;i++)
            dt[i]=(CFL*dx)/(maxUc);
        /*
	// Jameson Runge Kutta
        scalarF(Flux,F,U,c,W);
        //StegerWarmingF(Flux,W,U,c,rho);
        //correctedModifiedStegerWarmingF(Flux,W,U,c,rho,p);
        //modifiedStegerWarmingF(Flux,W,U,c,rho);
        //RoeF(Flux,W,F,U,c,rho,p);

        // Residual 0
        for(int k=0;k<3;k++)
        {
            for(int i=0;i<nx-1;i++)
            {
                Resi0[k][i+1]=Flux[k][i+1]*S[i+1]-Flux[k][i]*S[i]-Q[k][i+1];
            }
            Resi0[k][0]=0;
            Resi0[k][nx-1]=0;
        }
        // Runge Kutta 1
        for(int k=0;k<3;k++)
            for(int i=1;i<nx-1;i++)
            {
                W1[k][i]=W[k][i]-(dt[i]/2)*(Resi0[k][i])/dx;
            }

        for(int k=0;k<3;k++)
        {
            W1[k][0]=W[k][0];
            W1[k][nx-1]=W[k][nx-1];
        }
        for(int i=0;i<nx;i++)
        {
            Utemp[i]=W1[1][i]/W1[0][i];
            rhotemp[i]=W1[0][i];
            ptemp[i]=(gam-1)*(W1[2][i]-rhotemp[i]*Utemp[i]*Utemp[i]/2);
            ctemp[i]=sqrt((gam*ptemp[i]/rhotemp[i]));
            F1[0][i]=rhotemp[i]*Utemp[i];
            F1[1][i]=Utemp[i]*Utemp[i]*rhotemp[i]+ptemp[i];
            F1[2][i]=(W1[2][i]+ptemp[i])*Utemp[i];
            Q1[0][i]=0;
            Q1[2][i]=0;
        }
        for(int i=1;i<nx-1;i++)
            Q1[1][i]=ptemp[i]*(S[i]-S[i-1]);

        scalarF(Flux,F1,Utemp,ctemp,W1);
        //StegerWarmingF(Flux,W1,Utemp,ctemp,rhotemp);
        //correctedModifiedStegerWarmingF(Flux,W1,Utemp,ctemp,rhotemp,ptemp);
        //modifiedStegerWarmingF(Flux,W1,Utemp,ctemp,rhotemp);
        //RoeF(Flux,W1,F1,Utemp,ctemp,rhotemp,ptemp);

        // Residual 1
        for(int k=0;k<3;k++)
        {
            for(int i=0;i<nx-1;i++)
            {
                Resi1[k][i+1]=Flux[k][i+1]*S[i+1]-Flux[k][i]*S[i]-Q1[k][i+1];
            }
            Resi1[k][0]=0;
            Resi1[k][nx-1]=0;
        }
        // Step 2 Runge Kutta
        for(int k=0;k<3;k++)
            for(int i=1;i<nx-1;i++)
            {
                W2[k][i]=W[k][i]-(dt[i]/2)*(Resi1[k][i])/dx;
            }
        for(int k=0;k<3;k++)
        {
            W2[k][0]=W[k][0];
            W2[k][nx-1]=W[k][nx-1];
        }

        for(int i=0;i<nx;i++)
        {
            Utemp[i]=W2[1][i]/W2[0][i];
            rhotemp[i]=W2[0][i];
            ptemp[i]=(gam-1)*(W2[2][i]-rhotemp[i]*Utemp[i]*Utemp[i]/2);
            ctemp[i]=sqrt((gam*ptemp[i]/rhotemp[i]));
            F2[0][i]=rhotemp[i]*Utemp[i];
            F2[1][i]=Utemp[i]*Utemp[i]*rhotemp[i]+ptemp[i];
            F2[2][i]=(W2[2][i]+ptemp[i])*Utemp[i];
            Q2[0][i]=0;
            Q2[2][i]=0;
        }
        for(int i=1;i<nx-1;i++)
            Q2[1][i]=ptemp[i]*(S[i]-S[i-1]);

        scalarF(Flux,F2,Utemp,ctemp,W2);
        //StegerWarmingF(Flux,W2,Utemp,ctemp,rhotemp);
        //correctedModifiedStegerWarmingF(Flux,W2,Utemp,ctemp,rhotemp,ptemp);
        //modifiedStegerWarmingF(Flux,W2,Utemp,ctemp,rhotemp);
        //RoeF(Flux,W2,F2,Utemp,ctemp,rhotemp,ptemp);

       // Residual 2
        for(int k=0;k<3;k++)
        {
            for(int i=0;i<nx-1;i++)
            {
                Resi2[k][i+1]=Flux[k][i+1]*S[i+1]-Flux[k][i]*S[i]-Q2[k][i+1];
            }
            Resi2[k][0]=0;
            Resi2[k][nx-1]=0;
        }
        // Step 3 Runge Kutta
        for(int k=0;k<3;k++)
            for(int i=1;i<nx-1;i++)
            {
                W3[k][i]=W[k][i]-(dt[i]/2)*(Resi2[k][i])/dx;
            }

        for(int k=0;k<3;k++)
            for(int i=1;i<nx-1;i++)
            {
                W[k][i]=((double)1/6)*(W[k][i]+2*W1[k][i]+2*W2[k][i]+W3[k][i]);
            }
        for(int k=0;k<3;k++)
            for(int i=0;i<nx;i++)
            {
                Resi[k][i]=(2*Resi0[k][i]+2*Resi1[k][i]+Resi2[k][i])/6;
                //Resi[k][i]=(Resi2[k][i]);
                //printf("%f\n",Resi[k][i]);
            }
        */
        // Beggining of Euler explicit
        // Flux for each CV
        StegerWarmingF(Flux, W, U, c,rho);
        //modifiedStegerWarmingF(Flux, W, U, c,rho);
        //correctedModifiedStegerWarmingF(Flux,W,U,c,rho,p);
        //scalarF(Flux,F,U,c,W);
        //RoeF(Flux,W,F,U,c,rho,p);
//	for(int i=0;i<nx;i++)
//		printf("%f\n",Flux[1][i]);
        // Residual
        for(int k=0;k<3;k++)
        {
            for(int i=0;i<nx;i++)
            {
                Resi[k][i+1]=Flux[k][i+1]*S[i+1]-Flux[k][i]*S[i]-Q[k][i+1];
            }
            Resi[k][0]=0;
            Resi[k][nx-1]=0;
        }
        // Marching forward in time
        for(int k=0;k<3;k++)
            for(int i=1;i<nx-1;i++)
            {
                W[k][i]=W[k][i]-(dt[i]/V[i-1])*Resi[k][i];
            }
        // End of Euler explicit

        // Inlet Boundary Condition
        if(Mach[0]<1)
        {
            dpdu=pt*(gam/(gam-1))*pow(1-((gam-1)/(gam+1))*U[0]*U[0]/a2,1/(gam-1))*(-2*((gam-1)/(gam+1))*U[0]/a2);
            eigenvalues[0]=((U[1]+U[0]-c[1]-c[0])/2)*(dt[0]/dx);
            du=(-eigenvalues[0]*(p[1]-p[0]-rho[0]*c[0]*(U[1]-U[0])))/(dpdu-rho[0]*c[0]);

            U[0]=U[0]+du;
            T[0]=Tt*(1-((gam-1)/(gam+1))*U[0]*U[0]/a2);
            p[0]=pt*pow(T[0]/Tt,gam/(gam-1));
            rho[0]=p[0]/(R*T[0]);
            e[0]=rho[0]*(Cv*T[0]+0.5*U[0]*U[0]);
            c[0]=sqrt(gam*p[0]/rho[0]);
            Mach[0]=U[0]/c[0];

        }
        // Exit boundary condition
        // NOTE NOT SURE IF DX IS FULL DX OR DX/2
        eigenvalues[0]=((U[nx-1]+U[nx-2])/2)*(dt[nx-1]/(dx));
        eigenvalues[1]=((U[nx-1]+U[nx-2])/2+(c[nx-1]+c[nx-2])/2)*(dt[nx-1]/(dx));
        eigenvalues[2]=((U[nx-1]+U[nx-2])/2-(c[nx-1]+c[nx-2])/2)*(dt[nx-1]/(dx));
        charRel[0]=-eigenvalues[0]*(rho[nx-1]-rho[nx-2]-(1/(c[nx-1]*c[nx-1]))*(p[nx-1]-p[nx-2]));
        charRel[1]=-eigenvalues[1]*(p[nx-1]-p[nx-2]+rho[nx-1]*c[nx-1]*(U[nx-1]-U[nx-2]));
        charRel[2]=-eigenvalues[2]*(p[nx-1]-p[nx-2]-rho[nx-1]*c[nx-1]*(U[nx-1]-U[nx-2]));

        MachBound=((U[nx-1]+U[nx-2])/2)/((c[nx-1]+c[nx-2])/2);
        //MachBound=(U[nx-1]/(c[nx-1]));
        if(MachBound>1)
            dp=0.5*(charRel[1]+charRel[2]);
        else
            dp=0;

        drho=charRel[0]+dp/(pow(c[nx-1],2));
        du=(charRel[1]-dp)/(rho[nx-1]*c[nx-1]);

        U[nx-1]=U[nx-1]+du;
        rho[nx-1]=rho[nx-1]+drho;
        p[nx-1]=p[nx-1]+dp;
        T[nx-1]=p[nx-1]/(rho[nx-1]*R);
        e[nx-1]=rho[nx-1]*(Cv*T[nx-1]+0.5*pow(U[nx-1],2));
        c[nx-1]=sqrt((gam*p[nx-1])/rho[nx-1]);
        Mach[nx-1]=U[nx-1]/c[nx-1];

        // Update flow properties
        for(int i=1;i<nx-1;i++)
        {
            rho[i]=W[0][i];                                 // rho
            U[i]=W[1][i]/W[0][i];                           // U
            p[i]=(gam-1)*(W[2][i]-rho[i]*pow(U[i],2)/2);  // Pressure
            T[i]=p[i]/(rho[i]*R);                           // Temperature
            e[i]=W[2][i];                                   // Energy
            c[i]=sqrt((gam*p[i])/rho[i]);                 // Speed of sound
            Mach[i]=U[i]/c[i];                              // Mach number
        }
        // Update vectors
        for(int i=0;i<nx;i++)
        {
            W[0][i]=rho[i];
            F[0][i]=rho[i]*U[i];
            Q[0][i]=0;
            W[1][i]=rho[i]*U[i];
            F[1][i]=rho[i]*pow(U[i],2)+p[i];
            W[2][i]=e[i];
            F[2][i]=(e[i]+p[i])*U[i];
            Q[2][i]=0;
        }

        for(int i=1;i<nx-1;i++)
            Q[1][i]=p[i]*(S[i]-S[i-1]);

        // Calculating the norm of the density residual
        normm=0;
        for(int i=0;i<nx;i++)
            normm=normm+Resi[0][i]*Resi[0][i];
        normm=sqrt(normm);

        fprintf(DensityResiNorm,"%f,",normm);
        fprintf(DensityResiNorm,"\n");
        /*
        for(int i=0;i<nx;i++)
            fprintf(test,"%f,",Mach[i]);
        fprintf(test,"\n");
        /*
        for(int i=0;i<nx;i++)
            fprintf(test2,"%f,",Resi[0][i]);
        fprintf(test2,"\n");
        */
    }

    for(int i=0;i<nx;i++)
    {
            fprintf(PressureDist,"%f,",p[i]/pt);
        fprintf(PressureDist,"\n");
    }

    for(int i=0;i<nx;i++)
    {
            fprintf(EnergyDist,"%f,",W[2][i]);
        fprintf(EnergyDist,"\n");
    }

    for(int i=0;i<nx;i++)
    {
            fprintf(MachDist,"%f,",Mach[i]);
        fprintf(MachDist,"\n");
    }
	for(int k=0;k<3;k++)
	{
		std::cout<<"W"<<k+1<<std::endl;
		for(int i=0;i<nx;i++)
			std::cout<<W[k][i]<<std::endl;
	}



    fprintf(PressureLoss,"%f, %f",dx, (pt-((p[nx-1])/pow(1-(((gam-1)/(gam+1))*(pow(U[nx-1],2))/(a2*T[nx-1]/Tt)),(gam-1)/gam)))/pt);

    fclose(DensityResiNorm);
    fclose(PressureDist);
    fclose(MachDist);
    fclose(xData);
    fclose(IterationsData);
    fclose(PressureLoss);
    fclose(EnergyDist);
}
