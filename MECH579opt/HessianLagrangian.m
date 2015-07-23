function [ HL ] = HessianLagrangian( x,lambda )
%Hessian of the Lagrangian wrt to x
%HL=HessianRosen(x)-HessianCons(x);
HC= HessianCons(x);
HL=HessianRosen(x);
HL=HL-lambda(1)*HC(:,:,1);%-lambda(2)*HC(:,:,2);
end

