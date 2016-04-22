clc;clear;

n = 11;
nctl = n + 1;
deg = 2;
nknots = nctl + deg +1;
T = [zeros(deg,1);
     linspace(0, 1, nknots - 2 * deg)';
     ones(deg,1)];
u = linspace(0+1e-5,1-1e-5,100)';

for iu = 1:size(u,1)
    uy0(iu) = getbij(u(iu), 7, 2, T);   
    uy1(iu) = getbij(u(iu), 8, 2, T);
    uy2(iu) = getbij(u(iu), 9, 2, T);
    uy3(iu) = getbij(u(iu), 1, 2, T);
end

plot(u, uy0,'-o',u, uy1,'-o',u, uy2,'-o',u, uy3,'-o')
%%

x = linspace(1, 10, nctl)';
y = sin(x);

xx = zeros(100,1);
yy = zeros(100,1);

x(7) = x(7) + 2

k=2;
for iu = 1:size(u,1)
    for ib = 1 : nctl
        xx(iu) = xx(iu) + x(ib) * getbij(u(iu), ib, k, T);
        yy(iu) = yy(iu) + y(ib) * getbij(u(iu), ib, k, T);
    end
end
xx
yy
plot(xx, yy, '-o',x,y,'-o')
%%
clear; clc;

xstart = 0;
xend = 1;
xeps = 1e-15;

nx = 100;
x = linspace(0,1,nx)'
y = sin(10*x)
noise = 0.5 * rand(length(y), 1);
y = y + noise;
plot(x,y,'-o')


n = 9;
nctl = n + 1;
deg = 2;
nknots = nctl + deg +1;
T = [zeros(deg,1);
     linspace(0, 1+xeps, nknots - 2 * deg)';
     (xend + xeps) * ones(deg,1)];

for iu = 1:size(x,1)
    for ib = 1 : nctl
        A(iu, ib) = getbij(x(iu), ib, deg, T);
    end
end

ctls = A\y;

yy = zeros(nx,1);

for iu = 1:size(x,1)
    for ib = 1 : nctl
        yy(iu) = yy(iu) + ctls(ib) * getbij(x(iu), ib, deg, T);
    end
end

plot(x, yy, '-o',x,y,'o')
