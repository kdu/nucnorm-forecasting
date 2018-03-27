function w = froweights(L,K )
%froweights Number of times each element of an N-vector appears in an L x K
%Hankel matrix
N=L+K-1;
for i=1:L-1
    w(i)=i;
end

for i=L:K
    w(i)=L;
end

for i=K+1:N
    w(i)=N-i+1;
end

w=w';
end

