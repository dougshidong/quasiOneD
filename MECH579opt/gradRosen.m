function [yprime] = gradRosen(x)
%Gradient of the Rosenbrock function
yprime(1)=2*(200*x(1)^3-200*x(1)*x(2)+x(1)-1);
yprime(2)=200*(x(2)-x(1)^2);
yprime=yprime';
end

