#include<Eigen/Core>
#include<vector>
#include"globals.h"

Eigen::MatrixXd evaldSdDesign(
    std::vector <double> x,
    std::vector <double> dx,
    std::vector <double> designVar)
{
    Eigen::MatrixXd dSdDesign(nx + 1, designVar.size());
    if(desParam == 0) dSdDesign.setIdentity();
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
                dSdDesign(i, 0) = 0;
                dSdDesign(i, 1) = 0;
                dSdDesign(i, 2) = 0;
            }
            else
            {
                xh = fabs(x[i] - dx[i] / 2.0);
                dSdDesign(i, 0) = - pow(sin(PI * pow(xh, d2)), d3);
                dSdDesign(i, 1) = - d1 * d3 * PI * pow(xh, d2)
                                  * cos(PI * pow(xh, d2)) * log(xh)
                                  * pow(sin(PI * pow(xh, d2)), d3 - 1);
                dSdDesign(i, 2) = - d1 * log(sin(PI * pow(xh, d2)))
                                  * pow(sin(PI * pow(xh, d2)), d3);
            }
        }
    }
    return dSdDesign;
}
