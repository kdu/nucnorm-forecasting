load fort.mat;

Y=fort(1:110);
L=30;

X=hmat(Y,L);




val=0;
for j=1000:200:100000
    val=val+1;
    Ya(:,val)=mcwf(Y,L,10,froweights(L,size(X,2)),j);
end

N=length(Y)
for val=1:496
smse(val)=sqrt((1/10)*sum((Ya(N+1:N+10,val)-fort(N+1:N+10)).^2));
end