function [dfdx] = finiteD1(x)
h=0.000000001;
current=evalf(x);

dfdx(1)=(evalf([x(1)+h,x(2),x(3),x(4)])-current)/(h);
dfdx(2)=(evalf([x(1),x(2)+h,x(3),x(4)])-current)/(h);
dfdx(3)=(evalf([x(1),x(2),x(3)+h,x(4)])-current)/(h);
dfdx(4)=(evalf([x(1),x(2),x(3),x(4)+h])-current)/(h);
dfdx=dfdx';

end
