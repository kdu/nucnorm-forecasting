function [ y ] = hankvec( X )
%hankvec Read a vector from a Hankel matrix X

y=[X(1,1:end-1)';X(1:end,end)];

end

