function [c,A] = evalcons(x,y)
% Evaluate the constraint and its gradient at x
%Constraint
c(1)=x(3)-y(1);
c(2)=x(4)-y(2);
c=c';

%Gradient of the constraint
%First row, derivative of c1
A(1,1)=0;
A(1,2)=0;
A(1,3)=1;
A(1,4)=0;
%Second row, derivative of c2
A(2,1)=0;
A(2,2)=0;
A(2,3)=0;
A(2,4)=1;
%A=-A;
end

