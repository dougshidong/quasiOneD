function [dfdx] = finiteD2(x)
h=0.000000001;
dfdx(1)=(evalf([x(1)+h,x(2),x(3),x(4)])-evalf([x(1)-h,x(2),x(3),x(4)]))/(2*h);
dfdx(2)=(evalf([x(1),x(2)+h,x(3),x(4)])-evalf([x(1),x(2)-h,x(3),x(4)]))/(2*h);
dfdx(3)=(evalf([x(1),x(2),x(3)+h,x(4)])-evalf([x(1),x(2),x(3)-h,x(4)]))/(2*h);
dfdx(4)=(evalf([x(1),x(2),x(3),x(4)+h])-evalf([x(1),x(2),x(3),x(4)-h]))/(2*h);
dfdx=dfdx';

end

