function [ y ] = hankvec_avg( X )
%hankvec Perform diagonal averaging 
  [L,K] = size(X);
  Xvec = reshape([zeros(K);X], (K+L)*K, 1);
  Bmat = reshape(Xvec(K+1:end), (K+L-1), K);
  y = sum(Bmat,2)./froweights(L,K);
end

