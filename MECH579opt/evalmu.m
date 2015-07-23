function [ mmu ] = evalmu(grad,pk,HL,c,rho)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%rho=0.0000001;
sigma=0;
if(pk'*HL*pk>0) 
    sigma=0.5;
end

mmu=(grad'*pk+sigma*pk'*HL*pk)/((1-rho)*norm(c,1));
end

