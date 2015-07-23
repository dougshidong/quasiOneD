function [ L ] = Lagrangian( x,lambda)
%Lagrangian of the function
L=Rosenbrock(x)-lambda'*c1(x)
end

