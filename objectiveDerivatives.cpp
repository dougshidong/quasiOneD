#include<vector>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include"globals.h"
#include"quasiOneD.h"
#include<iostream>

using namespace Eigen;

// dIc / dW
VectorXd evaldIcdW(
    std::vector <double> W,
    std::vector <double> dx)
{
    VectorXd dIcdW(3 * nx);

    std::vector <double> ptarget(nx, 0);
    double dpdw[3], rho, u, p;
    ioTargetPressure(-1, ptarget);
    for(int i = 0; i < nx; i++)
    {
        rho = W[i * 3 + 0];
        u = W[i * 3 + 1] / rho;
        p = (gam - 1) * ( W[i * 3 + 2] - rho * u * u / 2.0 );

        dpdw[0] = (gam - 1) / 2.0 * u * u;
        dpdw[1] = - (gam - 1) * u;
        dpdw[2] = (gam - 1);

        dIcdW[i * 3 + 0] = (p / ptin - ptarget[i]) * dpdw[0] * dx[i] / ptin;
        dIcdW[i * 3 + 1] = (p / ptin - ptarget[i]) * dpdw[1] * dx[i] / ptin;
        dIcdW[i * 3 + 2] = (p / ptin - ptarget[i]) * dpdw[2] * dx[i] / ptin;
    }
    return dIcdW;
}

// ddIc / dWdW
SparseMatrix <double> evaldIcdWdW(
    std::vector <double> W,
    std::vector <double> dx)
{
    SparseMatrix <double> ddIcdWdW(3 * nx, 3 * nx);
    Matrix3d ddIcdWdW_temp = Matrix3d::Zero();
    Matrix3d ddpdWdW = Matrix3d::Zero();
    Vector3d dpdW = Vector3d::Zero();

    std::vector <double> ptarget(nx, 0);
    double rho, u, p;
    double dxptin2;
    ioTargetPressure(-1, ptarget);
    for(int Wi = 0; Wi < nx; Wi++)
    {
        rho = W[Wi * 3 + 0];
        u = W[Wi * 3 + 1] / rho;
        p = (gam - 1) * ( W[Wi * 3 + 2] - rho * u * u / 2.0 );

        dpdW(0) = (gam - 1) / 2.0 * u * u;
        dpdW(1) = - (gam - 1) * u;
        dpdW(2) = (gam - 1);
        
        ddpdWdW(0, 0) = -u * u * (gam - 1.0) / rho;
        ddpdWdW(1, 1) = (1.0 - gam) / rho;
        ddpdWdW(1, 0) = u * (gam - 1.0) / rho;
        ddpdWdW(0, 1) = ddpdWdW(1, 0);

        dxptin2 = dx[Wi] / pow(ptin, 2);
        ddIcdWdW_temp = dxptin2 * dpdW * dpdW.transpose()
                        + (dxptin2 * p - ptarget[Wi] * dx[Wi] / ptin) * ddpdWdW;

        for(int ki = 0; ki < 3; ki++)
        {
            for(int kj = 0; kj < 3; kj++)
            {
                ddIcdWdW.insert(Wi*3+ki, Wi*3+kj) = ddIcdWdW_temp(ki, kj);
            }
        }
    }
    return ddIcdWdW;
}

VectorXd evaldIcdS()
{
    VectorXd dIcdS(nx + 1);
    dIcdS.setZero();
    return dIcdS;
}

MatrixXd evalddIcdSdS()
{
    MatrixXd ddIcdSdS(nx + 1, nx + 1);
    ddIcdSdS.setZero();
    return ddIcdSdS;
}

MatrixXd evalddIcdWdS()
{
    MatrixXd ddIcdWdS(3 * nx, nx + 1);
    ddIcdWdS.setZero();
    return ddIcdWdS;
}
