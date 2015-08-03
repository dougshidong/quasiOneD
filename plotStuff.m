close all;

nx=200;

fid=fopen('Results.dat','r');
formatSpec='%f';

x=fscanf(fid,formatSpec,nx);
pressure=fscanf(fid,formatSpec,nx);
rho=fscanf(fid,formatSpec,nx);
Mach=fscanf(fid,formatSpec,nx);

xhalf=fscanf(fid,formatSpec,nx+1);
S=fscanf(fid,formatSpec,nx+1);

conv=fscanf(fid,formatSpec,inf);

figure(1)
plot(x,pressure,'-o');
title('Pressure Distribution');

figure(2)
plot(x,rho,'-o');
title('Density Distribution');

figure(3)
plot(x,Mach,'-o');
title('Mach Distribution');

figure(4)
plot(xhalf,S,'-o');
title('Geom');

figure(5)
loglog(conv(1:size(conv,1)/2), conv(size(conv,1)/2+1:size(conv,1)),'-+');
title('conv')

fclose(fid);