function [ sigma ] = sigmaFletcherR( r1,r2 )
%Calculation of sigma for with the Fletcher-Reeves method

sigma=(r2'*r2)/(r1'*r1);

end

