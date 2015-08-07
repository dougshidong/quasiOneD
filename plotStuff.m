close all;

fid=fopen('Results.dat','r');
formatSpec='%f';

nx=fscanf(fid,'%d',1);
x=fscanf(fid,formatSpec,nx);
pressure=fscanf(fid,formatSpec,nx);
rho=fscanf(fid,formatSpec,nx);
Mach=fscanf(fid,formatSpec,nx);

xhalf=fscanf(fid,formatSpec,nx+1);
S=fscanf(fid,formatSpec,nx+1);

conv=fscanf(fid,formatSpec,inf);

fclose(fid);

fid=fopen('targetP.dat','r');
nx=fscanf(fid,'%d',1);
targetp=fscanf(fid,formatSpec,nx);
fclose(fid);

figure
plot(x,pressure,'-o',x,targetp,'-x');
title('Pressure Distribution');
%%
figure
plot(x,rho,'-o');
title('Density Distribution');

figure
plot(x,Mach,'-o');
title('Mach Distribution');

figure
plot(xhalf,S,'-o');
title('Geom');

figure
loglog(conv(1:size(conv,1)/2), conv(size(conv,1)/2+1:size(conv,1)),'-+');
title('conv')

fclose(fid);