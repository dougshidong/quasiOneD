close all
clear;
global total
total=0;
method=4;
% 1 = STEEPEST
% 2 = NLCG (Fletcher-Reeves)
% 3 = NEWTON
% 4 = QUASI-NEWTON (BFGS)

diff=5;
% 1 = FINITE DIFFERENCE 1ST ORDER
% 2 = FINITE DIFFERENCE 2ND ORDER
% 3 = COMPLEX STEP
% 4 = DIRECT METHOD
% 5 = ADJOINT METHOD

n=4; % number of design variables
m=2; % number of equality constraints
i=1;
normz=0;

x1=[-0.1,-1,1,1]';
x2=x1(:,i);
lambda=[0,0];
if(m~=0)
    lambda(1:m)=0.25;
    lambda=lambda';
else
    lambda=[];
end
alpha=0.3;
tol=10^-12;

y=evaly(x1(:,i));

if(diff==1)
    gradz=(finiteD1(x1(:,i)));
elseif(diff==2)
    gradz=(finiteD2(x1(:,i)));
elseif(diff==3)
    gradz=(complexD(x1(:,i)));
else
    if(diff==4)
        gradz=directD(x1(:,i),y);
    elseif(diff==5)
        gradz=adjointD(x1(:,i),y);
    end
end

if(method==1 || method==2)
    pk=-gradz(:,i);
elseif(method==3)
    % Unconstrained
    %B=HessianRosen(x1(:,i));
    %pk=-inv(B)*gradz(:,i);
    
    % Constrained
    HL=HessianLagrangian(x1(:,i),lambda);
    [c,A]=evalcons(x1(:,i),y);
    Fprime=[HL -A';A zeros(m)];
    F=[gradz(:,i)-A'*lambda;c];
    pVector=-inv(Fprime)*F;
    pk=pVector(1:n);
    plambda=pVector(n+1:m+n);
    %alpha=stepBacktrack(x1(:,i),gradz(:,i),pk,HL,c); %unconstrained

elseif(method==4)
    % Unconstrained
    %B=eye(n)
    %H=eye(n);
    %pk=-H*gradz(:,i);
    %y=evaly(x1(:,i),yt);
    %alpha=stepBacktrackUncons(x1(:,i),y,pk);
    %pk=[1;1]
    %stepBacktrackUncons(x1(:,i),pk);    
    
    %Constrained
    B=eye(n);
    [c,A]=evalcons(x1(:,i),y);
    Fprime=[B -A';A zeros(m)];
    F=[gradz(:,i)-A'*lambda;c];
    pVector=-inv(Fprime)*F;
    pk=pVector(1:n);
    plambda=pVector(n+1:m+n);
    alpha=stepBacktrackCons(x1(:,i),gradz(:,i),pk,B,c);
    
end
normz=norm(gradz(:,i)-A'*lambda);
currentnorm=normz(i)

f=evalf(x1);

while(normz(i)>tol)    
    
    s=alpha*pk;
    x2=x1(:,i)+s;
    
    if(m~=0)
        lambda=lambda+alpha*plambda;
    end
    
    i=i+1;
    x1(:,i)=x2; 

    y=evaly(x1(:,i));

    if(diff==1)
        gradz(:,i)=(finiteD1(x1(:,i)));
    elseif(diff==2)
        gradz(:,i)=(finiteD2(x1(:,i)));
    elseif(diff==3)
        gradz(:,i)=(complexD(x1(:,i)));
    else
        if(diff==4)
            gradz(:,i)=directD(x1(:,i),y);
        elseif(diff==5)
            gradz(:,i)=adjointD(x1(:,i),y);
        end
    end

    if(method==1 || method==2)
        pk=-gradz(:,i);
    elseif(method==3)
        % Unconstrained
        %B=HessianRosen(x1(:,i));
        %pk=-inv(B)*gradz(:,i);

        % Constrained
        HL=HessianLagrangian(x1(:,i),lambda);
        [c,A]=evalcons(x1(:,i),y);
        Fprime=[HL -A';A zeros(m)];
        F=[gradz(:,i)-A'*lambda;c];
        pVector=-inv(Fprime)*F;
        pk=pVector(1:n);
        plambda=pVector(n+1:m+n);
        %alpha=stepBacktrack(x1(:,i),gradz(:,i),pk,HL,c); %unconstrained

    elseif(method==4)
        % Unconstrained
        %B=eye(n)
        %H=eye(n);
        %pk=-H*gradz(:,i);
        %y=evaly(x1(:,i),yt);
        %alpha=stepBacktrackUncons(x1(:,i),y,pk);
        %pk=[1;1]
        %stepBacktrackUncons(x1(:,i),pk);    

        %Constrained
        Atemp=A;        
        [c,A]=evalcons(x1(:,i),y);
        DG=(gradz(:,i)-A'*lambda)-(gradz(:,i-1)-Atemp'*lambda);
        rk=evalRK(s,B,DG,n);
        B=B+(rk*rk')/(s'*rk)-(B*s*s'*B)/(s'*B*s)
        %B=eye(n)
        Fprime=[B -A';A zeros(m)];
        F=[gradz(:,i)-A'*lambda;c];
        pVector=-inv(Fprime)*F;
        pk=pVector(1:n);
        plambda=pVector(n+1:m+n);
        alpha=stepBacktrackCons(x1(:,i),gradz(:,i),pk,B,c);
        if i>9
            alpha=1;
        end
    end
    normz(i)=norm(gradz(:,i)-A'*lambda);
    normz(i)=norm(gradz(1,i),gradz(2,i));
    currentnorm=normz(i);
    x2
    currentnorm
end
x2
iterations=1:i;

figure(1)
contourMatrix;
contour(CONTX,CONTY,CONTZ,75)
hold on
plot(x1(1,:),x1(2,:),'-o')
xlabel('x');
ylabel('y');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]); 

if(method==1)
    title('Steepest Descent Path')
elseif(method==2)
    title('Non-Linear Conjugate Gradient Path')
elseif(method==3)
    title('SQP Newton Path')
    %saveas(gcf, 'SQPNewtonPathc2_0m2', 'pdf')
elseif(method==4)
    title('Optimization Path')
    %saveas(gcf, 'PathIDF', 'pdf')
end

figure(2)
semilogy(iterations,normz)	
xlabel('Iterations');
ylabel('Gradient');
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
if(method==1)
    title('Steepest Descent Convergence')
elseif(method==2)
    title('Non-Linear Conjugate Gradient Convergence')
elseif(method==3)
    title('Gradient Convergence')
    %saveas(gcf, 'noisypath', 'pdf')
elseif(method==4)
    title('Gradient Convergence')
    %saveas(gcf, 'ConvIDF', 'pdf')
end