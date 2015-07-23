function [y] = evaly(x)
global total yi
y(2)=1;
y(1)=1;
change=1;
tol=10^-12;
yi=0;
while(change>tol)
    yi=yi+1;
    if(i>200)
        tol=10^-10;
    end
    if(i>500)
        tol=10^-5;
    end
    if(i>700)
        tol=1;
    end
    if(i>1000)
        tol=2;
    end
    oldy1=y(1);
    oldy2=y(2);
    y(1)=(-3*exp(-(x(1)+1)^2-.25*(x(2)+1)^2)+sin(y(2)));
    y(2)=(-3*exp(-5*(x(1)-3)^2-.25*(x(2)-3)^2)+exp(-y(1)));
    change=abs(y(1)-(-3*exp(-(x(1)+1)^2-.25*(x(2)+1)^2)+sin(y(2))));
    %sqrt(abs(oldy2-y(2))^2+abs(oldy1-y(1))^2);
end

total=total+yi;