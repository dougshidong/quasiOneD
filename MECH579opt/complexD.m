function [dfdx] = complexD(x)
%Complex step derivative
h(1)=10^-16;
dfdx(1)=imag(evalf([x(1)+1i*h,x(2),x(3),x(4)]))/h;
dfdx(2)=imag(evalf([x(1),x(2)+1i*h,x(3),x(4)]))/h;
dfdx(3)=imag(evalf([x(1),x(2),x(3)+1i*h,x(4)]))/h;
dfdx(4)=imag(evalf([x(1),x(2),x(3),x(4)+1i*h]))/h;
dfdx=dfdx';
end