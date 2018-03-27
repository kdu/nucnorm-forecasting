function [ X ] = hmat( Y, L )
%hmat Hankel matrix from a vector Y
%   The usual Hankel matrix construction from a vector in MATLAB might not
%   help if there are missing values in Y
N1=size(Y,1);
N2=size(Y,2);
N=max(N1,N2);
X=zeros(L,N+1-L);

for i=1:L
    for j=1:N+1-L
        X(i,j)=Y(i+j-1);
    end
end

end

