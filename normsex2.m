%compute the distance of the optimal rank r approximations
%of the deathsdata using alternative norms

N=72;
L=24;

load deathsdata.mat
Yfull=deathsdata;
Y=deathsdata(1:N);

X=hmat(Y,L);

X3=lra(X,3);
X6=lra(X,6);
X12=lra(X,12);

y3=[X3(1,1:size(X,2)),X3(2:end,end)']';
y6=[X6(1,1:size(X,2)),X6(2:end,end)']';
y12=[X12(1,1:size(X,2)),X12(2:end,end)']';

w=ones(N,1);
n31=norm(sqrt(w).*(y3-Y))
n61=norm(sqrt(w).*(y6-Y))
n121=norm(sqrt(w).*(y12-Y))

lambda=0.05;
for i=1:N
    w(i)=exp(lambda*i);
end

n31w=norm(sqrt(w).*(y3-Y))
n61w=norm(sqrt(w).*(y6-Y))
n121w=norm(sqrt(w).*(y12-Y))