function [ HC] = HessianCons(x)
%Hessian (wrt to x) of the constraint
HC(1:2,1:2,1)=zeros(2);
%HC(1:2,1:2,1)=-2*eye(2);
end

