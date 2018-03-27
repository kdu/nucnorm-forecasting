function [p] = sum_1d_exp( lambdas, len, coefs )
%SUM_1D_EXP Sum of 1D exponentials
%   Detailed explanation goes here
  r = length(lambdas);
  if nargin < 3 
    coefs = ones(r,1);
  end    
  p = ((ones(len,1) * lambdas) .^ ((0:len-1)' * ones(1,r))) * coefs;
end

