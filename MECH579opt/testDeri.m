clear;
x=[-0.1,-1,1,1];
y=evaly(x);
dm=directD(x,y);
am=adjointD(x,y);
for i=1:16
    h=10^-i;
    fd1(i)=(evalf([x(1)+h,x(2),x(3),x(4)])-evalf(x))/(h);
    fd2(i)=(evalf([x(1)+h,x(2),x(3),x(4)])-evalf([x(1)-h,x(2),x(3),x(4)]))/(2*h);
    cs(i)=imag(evalf([x(1)+1i*h,x(2),x(3),x(4)]))/h;
    %fd1(i)=abs(fd1(i)-am(1));
    %fd2(i)=abs(fd2(i)-am(1));
    %cs(i)=abs(cs(i)-am(1));
    hp(i)=10^-i;
end
fd1=fd1';
fd2=fd2';
cs=cs';

figure(1)
loglog(hp,fd1,'-o',hp,fd2,'-o',hp,cs,'-o')
set(gca,'XDir','reverse');
title('Derivative Comparison')
legend('FD1','FD2','CS')
xlabel('h');
ylabel('Error');
set(gcf, 'PaperPosition', [0 0 8 5]);
set(gcf, 'PaperSize', [8 5]);
%saveas(gcf, 'derivatives', 'pdf')