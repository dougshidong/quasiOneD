// Flux Splitting Schemes

#include<vector>
#include<math.h>
#include<iostream>
#include"globals.h"
#include"flux.h"
#include"convert.h"

// Get Flux based on FluxScheme
void getFlux(std::vector <double> &Flux,
             std::vector <double> W)
{
    if(FluxScheme == 0) // SW
        Flux_StegerWarming(Flux, W);
    else if(FluxScheme == 1) // Scalar
        Flux_Scalar(Flux, W);
    else if(FluxScheme == 2) // MWS
        Flux_MSW(Flux, W);
    else if(FluxScheme == 3) // CMWS
        Flux_CMSW(Flux, W);
    else if(FluxScheme == -1) // SW New?
        Flux_SW(Flux, W);
}
// Flux Jacobian for SW/MSW/CMSW
void Flux_Jacobian(
    std::vector <double> &Ap_list,
    std::vector <double> &An_list,
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

        rho = W[i * 3 + 0];
        u = W[i * 3 + 1] / rho;
        e = W[i * 3 + 2];
        c = sqrt( gam / rho * (gam - 1) * ( e - rho * u * u / 2 ) );
        S[0][0] = 1.0;
        S[1][0] = -u / rho;
        S[2][0] = 0.5 * u * u * beta;
        S[1][1] = 1.0 / rho;
        S[2][1] = -u * beta;
        S[2][2] = beta;
        Sinv[0][0] = 1.0;
        Sinv[1][0] = u;
        Sinv[2][0] = 0.5 * u * u;
        Sinv[1][1] = rho;
        Sinv[2][1] = u * rho;
        Sinv[2][2] = 1.0 / beta;
        C[0][0] = 1.0;
        C[1][1] = rho * c;
        C[2][1] = -rho * c;
        C[0][2] = -1.0 / (c * c);
        C[1][2] = 1.0;
        C[2][2] = 1.0;
        Cinv[0][0] = 1.0;
        Cinv[0][1] = 1.0 / (2 * c * c);
        Cinv[0][2] = 1.0 / (2 * c * c);
        Cinv[1][1] = 1.0 / (2 * rho * c);
        Cinv[1][2] = -1.0 / (2 * rho * c);
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
}

// StegerWarming
void Flux_StegerWarming(
    std::vector <double> &Flux,
    std::vector <double> W)
{
    std::vector <double> Ap_list(nx * 3 * 3, 0), An_list(nx * 3 * 3, 0);

    Flux_Jacobian(Ap_list, An_list, W);

    for(int i = 1; i < nx; i++)
    {
        Flux[i * 3 + 0] = 0;
        Flux[i * 3 + 1] = 0;
        Flux[i * 3 + 2] = 0;
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            int Ap_pos = ((i - 1) * 3 * 3) + (row * 3) + col;
            int An_pos = (i * 3 * 3) + (row * 3) + col;
            Flux[i * 3 + row] += Ap_list[Ap_pos] * W[(i - 1) * 3 + col]
                               + An_list[An_pos] * W[i * 3 + col];
        }
    }
}

// Modified StegerWarming
void Flux_MSW(
    std::vector <double> &Flux,
    std::vector <double> W)
{
    std::vector <double> Ap_list(nx * 3 * 3, 0), An_list(nx * 3 * 3, 0);

    std::vector <double> Wavg(3 * nx);

    for(int i = 0; i < nx - 1; i++)
    {
        for(int k = 0; k < 3; k++)
        {
            Wavg[i * 3 + k] = (W[i * 3 + k] + W[(i+1) * 3 + k]) / 2.0;
        }
    }
    Wavg[(nx-1)*3 + 0] = 1;
    Wavg[(nx-1)*3 + 1] = 1;
    Wavg[(nx-1)*3 + 2] = 1;


    Flux_Jacobian(Ap_list, An_list, Wavg);

    for(int i = 1; i < nx; i++)
    {
        Flux[i * 3 + 0] = 0;
        Flux[i * 3 + 1] = 0;
        Flux[i * 3 + 2] = 0;
        for(int row = 0; row < 3; row++)
        for(int col = 0; col < 3; col++)
        {
            int Ahalf_pos = ((i-1) * 3 * 3) + (row * 3) + col;
            Flux[i * 3 + row] += Ap_list[Ahalf_pos] * W[(i - 1) * 3 + col]
                               + An_list[Ahalf_pos] * W[i * 3 + col];
        }
    }
}

void Flux_CMSW(
    std::vector <double> &Flux,
    std::vector <double> W)
{
    std::vector <double> FluxSW(3 * (nx + 1)), FluxMSW(3 * (nx + 1));
    std::vector <double> Whalf(3 * nx);
    std::vector <double> p(nx);

    Flux_StegerWarming(FluxSW, W);
    Flux_MSW(FluxMSW, W);

    for(int i = 0; i < nx; i++)
    {
        p[i] = (gam - 1.0) * ( W[i * 3 + 2] - (pow(W[i * 3 + 1], 2.0) / W[i * 3 + 0]) / 2.0 );
    }

    for(int i = 0; i < nx - 1; i++)
    {
        for(int k = 0; k < 3; k++)
        {
            Whalf[i * 3 + k]= 1.0 / (1.0 + pow( (p[i+1] - p[i]) / std::min(p[i+1], p[i]), 2 ));
        }
    }

    for(int i = 1; i < nx; i++)
    {
        for(int k = 0; k < 3; k++)
        {
            Flux[i * 3 + k] = 0;
            Flux[i * 3 + k] = Whalf[(i - 1) * 3 + k] * FluxMSW[i * 3 + k]
                     + (1.0 - Whalf[(i - 1) * 3 + k]) * FluxSW[i * 3 + k];
        }
    }
}

std::vector <double> F(3 * nx);
void Flux_Scalar(
    std::vector <double> &Flux,
    std::vector <double> W)
{
    int ki, kim;
    double lambda;
    double avgu, avgc;
    double rho, e;
    std::vector <double> u(nx), c(nx);

    // Get Convective Variables
    WtoF(W, F);

    for(int i = 0; i < nx; i++)
    {
        rho = W[i * 3 + 0];
        e = W[i * 3 + 2];
        u[i] = W[i * 3 + 1] / rho;
        c[i] = sqrt( gam / rho * (gam - 1) * ( e - rho * u[i] * u[i] / 2.0 ) );
    }
    for(int i = 1; i < nx; i++)
    {
        avgu = ( u[i - 1] + u[i] ) / 2.0;
        avgc = ( c[i - 1] + c[i] ) / 2.0;
        lambda = std::max( std::max( fabs(avgu), fabs(avgu + avgc) ),
                           fabs(avgu - avgc) );
        lambda = avgu + avgc;

        for(int k = 0; k < 3; k++)
        {
            ki = i * 3 + k;
            kim = (i - 1) * 3 + k;
            Flux[ki] = 0.5 * (F[kim] + F[ki]) - 0.5 * Scalareps * lambda * (W[ki] - W[kim]);
        }
    }
}

void partialJacobian(
    double *pA,
    double rho, double u, double c,
    double lamb1, double lamb2, double lamb3)
{
    double beta = gam - 1.0;
    double uu = u * u;
    double cc = c * c;

    pA[0] = (4.0 * cc * lamb1 + 2.0 * c * (-lamb2 + lamb3) * u
            + beta * (-2.0 * lamb1 + lamb2 + lamb3) * uu)/(4.0 * cc);
    pA[1] = (c * (lamb2 - lamb3)
            + beta * (2.0 * lamb1 - lamb2 - lamb3) * u) / (2.0 * cc);
    pA[2] = (beta * (-2.0 * lamb1 + lamb2 + lamb3))/(2.0 * cc);
    pA[3] = (u * (cc * (4.0 * lamb1 - 2.0 * (lamb2 + lamb3)) +
            (-2.0 + beta) * c * (lamb2 - lamb3) * u +
            beta * (-2.0 * lamb1 + lamb2 + lamb3) * uu)) / (4.0 * cc);
    pA[4] = (cc * (lamb2 + lamb3) - (-1 + beta) * c * (lamb2 - lamb3) * u
            + beta * (2.0 * lamb1 - lamb2 - lamb3) * uu)/(2.0 * cc);
    pA[5] = (beta * (c * (lamb2 - lamb3)
            + (-2.0 * lamb1 + lamb2 + lamb3) * u)) / (2.0 * cc);
    pA[6] = (u * (-4 * cc * c * (lamb2 - lamb3)
            + 2.0 * beta * cc * (2.0 * lamb1 - lamb2 - lamb3) * u
            + 2.0 * (-1.0 + beta) * beta * c * (lamb2 - lamb3) * uu
            + beta * beta * (-2.0 * lamb1 + lamb2 + lamb3) * uu * u))
            / (8.0 * beta * cc);
    pA[7] = (2.0 * cc * c * (lamb2 - lamb3)
            - beta * (-1.0 + 2.0 * beta) * c * (lamb2 - lamb3) * uu
            + beta * beta * (2.0 * lamb1 - lamb2 - lamb3) * uu * u)
            / (4.0 * beta * cc);
    pA[8] = (2.0 * cc * (lamb2 + lamb3)
            + 2.0 * beta * c * (lamb2 - lamb3) * u
            + beta * (-2.0 * lamb1 + lamb2 + lamb3) * uu)
            / (4.0 * cc);
}

void Flux_SW(
    std::vector <double> &Flux,
    std::vector <double> W)
{
    double Ap_list[nx][9], An_list[nx][9];

    double rho, u, e, c;
    double eig[3], eigp[3], eign[3];

    for(int i = 0; i < nx; i++)
    {
        rho = W[i * 3 + 0];
        u = W[i * 3 + 1] / rho;
        e = W[i * 3 + 2];
        c = sqrt( gam / rho * (gam - 1) * ( e - rho * u * u / 2 ) );

        double eps = 0.1;
        eig[0] = u;
        eig[1] = u + c;
        eig[2] = u - c;
        for(int k = 0; k < 3; k++)
        {
            if(eig[k] > 0) eigp[k] = (eig[k] + sqrt(pow(eig[k], 2) + pow(eps, 2))) / 2.0;
            else eign[k] = (eig[k] - sqrt(pow(eig[k], 2) + pow(eps, 2))) / 2.0;
        }

        partialJacobian(Ap_list[i], rho, u, c, eigp[0], eigp[1], eigp[2]);
        partialJacobian(An_list[i], rho, u, c, eign[0], eign[1], eign[2]);
    }


    for(int i = 1; i < nx; i++)
    {
        Flux[i * 3 + 0] = 0;
        Flux[i * 3 + 1] = 0;
        Flux[i * 3 + 2] = 0;
        for(int row = 0; row < 3; row++)
        {
            for(int col = 0; col < 3; col++)
            {
                Flux[i * 3 + row] += Ap_list[i-1][row * 3 + col] * W[(i - 1) * 3 + col]
                                   + An_list[i][row * 3 + col] * W[i * 3 + col];
            }
        }
    }
}

