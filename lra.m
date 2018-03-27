function [Xout]=lra(X,r)
    [U,S,V]=svd(X,'econ');
    for i=r+1:size(S,1)
        S(i,i)=0;
    end
    Xout=U*S*V';
end