function [ rk ] = evalRK(s,B,DG,n)
%Evaluates rk for the constrained BFGS
if(s'*DG<0.2*s'*B*s)
    theta=0.8*(s'*B*s)/(s'*B*s-s'*DG);
else
    theta=1;
end
rk=theta*DG-(1-theta)*B*s;

end

