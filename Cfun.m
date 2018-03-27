function y = Cfun(rho,n,L,K)
%Cfun needs to be smaller than one for the nuclear norm
%   heuristic give a solution which co-incides with the rank minimization
%   problem for the rank-one case
%
%   Cfun(rho,n,L,K) 

m=L+K-1-n;

if abs(rho) ~= 1
    y=abs(rho^(2*(m+1)-1)*(rho^n))/sqrt(abs((rho^(2*L)-1)*(rho^(2*K)-1)));
end

if abs(rho)==1
    y=(m+1)/sqrt(L*K);
end

