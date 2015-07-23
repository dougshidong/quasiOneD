function [alpha] = stepBacktrackUncons(x,y,p)
%Calculate the step length by backtracking
c=10^-4;
alpha=1;
%Wolfe1=Rosenbrock(x)+c*alpha*gradRosen(x)'*p;
while(evalf(x+alpha*p)>(evalf(x)+c*alpha*directD(x,y)'*p))
    alpha=alpha*0.9;
end

