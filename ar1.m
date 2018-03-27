function [ y ] = ar1( n, rho)
%AR1 Generate ar(1) noise with unit variance 
    x = sqrt(1-rho^2)*randn(n,1);
    y = filter(1,[1 -rho],x);
end

