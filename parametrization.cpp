#include<Eigen/Core>
#include<vector>
#include<iostream>
#include"globals.h"
#include"grid.h"

using namespace Eigen;

MatrixXd evaldSdDes(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> designVar)
{
    MatrixXd dSdDes(nx + 1, nDesVar);
    if(desParam == 0) 
    {
        dSdDes.setZero();
        for(int Si = 1; Si < nx; Si++)
        {
            dSdDes(Si, Si - 1) = 1;
        }
    }
    if(desParam == 1)
    {
        double d1 = designVar[0];
        double d2 = designVar[1];
        double d3 = designVar[2];
        double xh;
        for(int i = 0; i < nx + 1; i++)
        {
            if(i == 0 || i == nx)
            {
                dSdDes(i, 0) = 0;
                dSdDes(i, 1) = 0;
                dSdDes(i, 2) = 0;
            }
            else
            {
                xh = fabs(x[i] - dx[i] / 2.0);
                dSdDes(i, 0) = - pow(sin(PI * pow(xh, d2)), d3);
                dSdDes(i, 1) = - d1 * d3 * PI * pow(xh, d2)
                                  * cos(PI * pow(xh, d2)) * log(xh)
                                  * pow(sin(PI * pow(xh, d2)), d3 - 1);
                dSdDes(i, 2) = - d1 * log(sin(PI * pow(xh, d2)))
                                  * pow(sin(PI * pow(xh, d2)), d3);
            }
        }
    }
    return dSdDes;
}

MatrixXd evalddSdDesdDes_FD(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> designVar,
    int Si)
{
    double pertdi, pertdj;
    double h = 1e-3;
    std::vector <double> S(nx + 1);
    std::vector <double> S1(nx+1), S2(nx+1), S3(nx+1), S4(nx+1);
    std::vector <double> designVard(3);
    S = evalS(designVar, x, dx, 1);

    MatrixXd ddSdDesdDes(3,3);
    for(int di = 0; di < 3; di++)
    {
        for(int dj = 0; dj < 3; dj++)
        {
            pertdi = designVar[di] * h;
            pertdj = designVar[dj] * h;

            for(int m = 0; m<3; m++)
            {
                designVard[m] = designVar[m];
            }

            if(di != dj)
            {
                // ddFluxdWdW1 = ddFlux(i)/dW(i-1)dW(i-1)
                // R1
                designVard[di] = designVar[di] + pertdi;
                designVard[dj] = designVar[dj] + pertdj;

                S1 = evalS(designVard, x, dx, 1);

                // R2
                designVard[di] = designVar[di] + pertdi;
                designVard[dj] = designVar[dj] - pertdj;

                S2 = evalS(designVard, x, dx, 1);

                // R3
                designVard[di] = designVar[di] - pertdi;
                designVard[dj] = designVar[dj] + pertdj;

                S3 = evalS(designVard, x, dx, 1);

                // R4
                designVard[di] = designVar[di] - pertdi;
                designVard[dj] = designVar[dj] - pertdj;

                S4 = evalS(designVard, x, dx, 1);


                ddSdDesdDes(di, dj) = 
                    (S1[Si] - S2[Si] - S3[Si] + S4[Si]) / (4.0 * pertdi * pertdj);

                designVard[di] = designVar[di];
                designVard[dj] = designVar[dj];
            }
            else
            {
                // R1
                designVard[di] = designVar[di] + 2.0 * pertdi;

                S1 = evalS(designVard, x, dx, 1);

                // R2
                designVard[di] = designVar[di] + pertdi;

                S2 = evalS(designVard, x, dx, 1);
                // R3
                designVard[di] = designVar[di] - pertdi;

                S3 = evalS(designVard, x, dx, 1);
                // R4
                designVard[di] = designVar[di] - 2.0 * pertdi;
                S4 = evalS(designVard, x, dx, 1);

                ddSdDesdDes(di, dj) =
                                (-S1[Si] + 16*S2[Si] - 30*S[Si] + 16*S3[Si] - S4[Si])
                                             / (12 * pertdi * pertdj);
            }
        }
    }
    return ddSdDesdDes;

}
std::vector <MatrixXd> evalddSdDesdDes(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> designVar)
{
    int nD = designVar.size();
    std::vector <MatrixXd> ddSdDesdDes(nx + 1);
    MatrixXd ddSidDesdDes(nD, nD);
    if(desParam == 0)
    {
        ddSidDesdDes.setZero();
        for(int Si = 0; Si < nx + 1; Si++)
        {
            ddSdDesdDes[Si] = ddSidDesdDes;
        }
    }
    if(desParam == 1)
    {
        double d1 = designVar[0];
        double d2 = designVar[1];
        double d3 = designVar[2];
        double xh;
        double xd2, spxd2, cpxd2 ;
        for(int Si = 0; Si < nx + 1; Si++)
        {
            if(Si == 0 || Si == nx)
            {
                ddSdDesdDes[Si] = Matrix3d::Zero();
            }
            else
            {
                xh = fabs(x[Si] - dx[Si] / 2.0);
                xd2 = pow(xh, d2);
                spxd2 = sin(PI * xd2);
                cpxd2 = cos(PI * xd2);
                // Diagonal
                ddSidDesdDes(0, 0) = 0.0;
                ddSidDesdDes(1, 1) = -0.5 * d1 * d3 * PI * xd2 * pow(log(xh), 2.0)
                                    * pow(spxd2, d3-2.0)
                                    * ((d3-2.0) * PI * xd2 
                                       + d3 * PI * xd2  * cos(2.0 * PI * xd2)
                                       + sin(2.0 * PI * xd2));
                ddSidDesdDes(2, 2) = -d1 * pow(log(spxd2), 2.0) * pow(spxd2, d3);

                // Off-Diagonal
                ddSidDesdDes(0, 1) = -d3 * PI * xd2 * cpxd2 * log(xh) * pow(spxd2, d3-1.0);
                ddSidDesdDes(1, 0) = ddSidDesdDes(0, 1);

                ddSidDesdDes(0, 2) = -log(spxd2) * pow(spxd2, d3);
                ddSidDesdDes(2, 0) = ddSidDesdDes(0, 2);

                ddSidDesdDes(1, 2) = -d1 * PI * xd2 * cpxd2 * log(xh)
                                    * (1.0 + d3 * log(spxd2))
                                    * pow(spxd2, d3-1.0);
                ddSidDesdDes(2, 1) = ddSidDesdDes(1, 2);
                
                ddSdDesdDes[Si] = ddSidDesdDes;
            }

            //std::cout<<"Si "<<Si<<std::endl;
            //std::cout<<(ddSidDesdDes 
            //    - evalddSdDesdDes_FD(x, dx, designVar, Si)).norm()<<std::endl;
        }
    }
    return ddSdDesdDes;
}

