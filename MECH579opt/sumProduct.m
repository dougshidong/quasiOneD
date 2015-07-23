function [lambdaGradC] = sumProduct(lambda,A)
% Sums the product of lamda(i) and the gradC(i) from 1 to m
lambdaGradC=0;
for i=1:size(A,1)
    lambdaGradC=lambdaGradC+lambda(i)*A(1,:)'+lambda(i)*A(2,:)';
end
lambdaGradC=A'*lambda;
end
