function [alpha] = stepBacktrackCons(x,grad,pk,HL,c)
%Calculate the step length by backtracking
global total
rho=10^-4;
alpha=1;
mmu=evalmu(grad,pk,HL,c,rho);
muc=mmu*norm(c,1);
current=evalf(x);
total=total-2;
while((evalf(x+alpha*pk)+muc)>(current+muc)+alpha*rho*(-pk'*grad-muc) && alpha>10^-2)
    alpha=alpha*0.8;
    total=total-1;
end

end