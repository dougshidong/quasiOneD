function [y] = evaly(x)
global total;

y(1)=(-3*exp(-(x(1)+1)^2-.25*(x(2)+1)^2)+sin(x(4)));
y(2)=(-3*exp(-5*(x(1)-3)^2-.25*(x(2)-3)^2)+exp(-x(3)));

total=total+1;

end