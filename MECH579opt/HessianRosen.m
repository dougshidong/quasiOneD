function [ Hessian ] = HessianRosen(x)
%Calculates the Hessian of the Rosenbrock function
Hessian(1,1)=1200*x(1)^2-400*x(2)+2;
Hessian(2,2)=200;
Hessian(1,2)=-400*x(1);
Hessian(2,1)=-400*x(1);
end

