// Flux Splitting Schemes

#include<vector>
#include<math.h>
#include<iostream>
#include"globals.h"
#include"flovar.h"
#include"flux.h"


// Get Flux based on FluxScheme
void getFlux(std::vector <double> &Flux,
             std::vector <double> W,
             std::vector <double> F)
{
    if(FluxScheme == 0) // SW
        Flux_StegerWarming(Flux, W);
    else if(FluxScheme == 1) // Scalar
        Flux_Scalar(Flux, W, F);
}

// StegerWarming
void Flux_StegerWarming(std::vector <double> &Flux,
                         std::vector <double> W)
{
    double eps = 0.1;

    double S[3][3] = {{0}},
           Sinv[3][3] = {{0}},
           C[3][3] = {{0}},
           Cinv[3][3] = {{0}},
           lambdaP[3][3],
           lambdaN[3][3];
    double lambdaa[3];
    
    
    double Ap[3][3], An[3][3], tempP[3][3], tempN[3][3], prefix[3][3], suffix[3][3];

    std::vector <double> Ap_list(nx * 3 * 3, 0), An_list(nx * 3 * 3, 0);

    double beta = 0.4;//gam-1;
    double rho, u, e, c;

    for(int i = 0; i < nx; i++)
    {
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            Ap[row][col] = 0;
            An[row][col] = 0;
            tempP[row][col] = 0;
            tempN[row][col] = 0;
            prefix[row][col] = 0;
            suffix[row][col] = 0;
            lambdaP[row][col] = 0;
            lambdaN[row][col] = 0;
        }
        
        rho = W[0 * nx + i];
        u = W[1 * nx + i] / rho;
        e = W[2 * nx + i];
        c = sqrt( gam / rho * (gam - 1) * ( e - rho * u * u / 2 ) );
        S[0][0] = 1;
        S[1][0] = -u / rho;
        S[2][0] = 0.5 * u * u * beta;
        S[1][1] = 1 / rho;
        S[2][1] = -u * beta;
        S[2][2] = beta;
        Sinv[0][0] = 1;
        Sinv[1][0] = u;
        Sinv[2][0] = 0.5 * u * u;
        Sinv[1][1] = rho;
        Sinv[2][1] = u * rho;
        Sinv[2][2] = 1 / beta;
        C[0][0] = 1;
        C[1][1] = rho * c;
        C[2][1] = -rho * c;
        C[0][2] = -1 / (c * c);
        C[1][2] = 1;
        C[2][2] = 1;
        Cinv[0][0] = 1;
        Cinv[0][1] = 1 / (2 * c * c);
        Cinv[0][2] = 1 / (2 * c * c);
        Cinv[1][1] = 1 / (2 * rho * c);
        Cinv[1][2] = -1 / (2 * rho * c);
        Cinv[2][1] = 0.5;
        Cinv[2][2] = 0.5;
        lambdaa[0] = u;
        lambdaa[1] = u + c;
        lambdaa[2] = u - c;
        
        for(int k = 0; k < 3; k++)
            if(lambdaa[k] > 0)
                lambdaP[k][k] = (lambdaa[k] + sqrt(pow(lambdaa[k], 2) + pow(eps, 2))) /2;
            else
                lambdaN[k][k] = (lambdaa[k] - sqrt(pow(lambdaa[k], 2) + pow(eps, 2))) / 2;

        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                prefix[row][col] += Sinv[row][k] * Cinv[k][col];
                suffix[row][col] += C[row][k] * S[k][col];
            }
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                tempP[row][col] += prefix[row][k] * lambdaP[k][col];
                tempN[row][col] += prefix[row][k] * lambdaN[k][col];
            }
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
            for(int k = 0; k < 3; k++)
            {
                Ap[row][col] += tempP[row][k] * suffix[k][col];
                An[row][col] += tempN[row][k] * suffix[k][col];
            }
        // could remove above loop and just use aplist and anlist
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            int vec_pos = (i * 3 * 3) + (row * 3) + col;
            Ap_list[vec_pos] = Ap[row][col];
            An_list[vec_pos] = An[row][col];
        }

    }

    for(int i = 1; i < nx; i++)
    {
        Flux[0 * nx + i] = 0;
        Flux[1 * nx + i] = 0;
        Flux[2 * nx + i] = 0;
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            int Ap_pos = ((i - 1) * 3 * 3) + (row * 3) + col;
            int An_pos = (i * 3 * 3) + (row * 3) + col;
            Flux[row * nx + i] = Flux[row * nx + i] + Ap_list[Ap_pos] * W[col * nx + (i - 1)]
                 + An_list[An_pos] * W[col * nx + i];
        }
    }

}


void Flux_Scalar(std::vector <double> &Flux,
                 std::vector <double> W,
                 std::vector <double> F)
{
    int ki;
    double eps = 0.5, lambda;
    double avgu, avgc;
    double rho, e;
    std::vector <double> u(nx), c(nx);

    for(int i = 0; i < nx; i++)
    {
        rho = W[0 * nx + i];
        e = W[2 * nx + i];
        u[i] = W[1 * nx + i] / rho;
        c[i] = sqrt( gam / rho * (gam - 1) * ( e - rho * u[i] * u[i] / 2 ) );
    }
    for(int i = 1; i < nx; i++)
    {
        avgu = ( u[i - 1] + u[i] ) / 2;
        avgc = ( c[i - 1] + c[i] ) / 2;
        lambda = std::max( std::max( fabs(avgu), fabs(avgu + avgc) ),
                           fabs(avgu - avgc) );

        for(int k = 0; k < 3; k++)
        {
            ki = k * nx + i;
            Flux[ki] = 0.5 * (F[ki - 1] + F[ki]) - 0.5 * eps * lambda * (W[ki] - W[ki - 1]);
        }
    }
}

