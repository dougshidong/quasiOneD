function [dfdx] = directD(x,y)

%Partial derivatives of f
pfpx(1)=-(20*(-2*x(1)+2))*exp(-(x(1)-1)^2-.25*(x(2)-1)^2);
pfpx(2)=-(20*(-.50*x(2)+.50))*exp(-(x(1)-1)^2-.25*(x(2)-1)^2);
pfpx(3)=0;
pfpx(4)=0;
pfpy1=1;
pfpy2=-sin(y(2));

%Partial derivatives of R1
pR1px(1)=-(3*(-2*x(1)-2))*exp(-(x(1)+1)^2-.25*(x(2)+1)^2);
pR1px(2)=-(3*(-.50*x(2)-.50))*exp(-(x(1)+1)^2-.25*(x(2)+1)^2);
pR1px(3)=0;
pR1px(4)=cos(x(4));
pR1py1=-1;
pR1py2=0;

%Partial derivatives of R2
pR2px(1)=-(3*(-10*x(1)+30))*exp(-5*(x(1)-3)^2-.25*(x(2)-3)^2);
pR2px(2)=-(3*(-.50*x(2)+1.50))*exp(-5*(x(1)-3)^2-.25*(x(2)-3)^2);
pR2px(3)=-exp(-x(3));
pR2px(4)=0;
pR2py1=0;
pR2py2=-1;

%Form and inverse matrix A
A=[pR1py1 pR1py2; pR2py1 pR2py2];
Ainv=inv(A);

%Calculate dydx
dydx(:,1)=-Ainv*[pR1px(1);pR2px(1)];
dydx(:,2)=-Ainv*[pR1px(2);pR2px(2)];
dydx(:,3)=-Ainv*[pR1px(3);pR2px(3)];
dydx(:,4)=-Ainv*[pR1px(4);pR2px(4)];

%Calculate dfdx using pfp* and dydx(*,*)
dfdx(1)=pfpx(1)+pfpy1*dydx(1,1)+pfpy2*dydx(2,1);
dfdx(2)=pfpx(2)+pfpy1*dydx(1,2)+pfpy2*dydx(2,2);
dfdx(3)=pfpx(3)+pfpy1*dydx(1,3)+pfpy2*dydx(2,3);
dfdx(4)=pfpx(4)+pfpy1*dydx(1,4)+pfpy2*dydx(2,4);

dfdx=dfdx';
end