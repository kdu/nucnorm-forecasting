function [ d] = wvnorm( y,w )
%wvnorm Compute weighted vector norm given a vector y and a 
% vector of weights w

d=sum(w.*(y.^2));

end

