load fort.mat;

Y=fort;
L=60;

X=hmat(fort,L);
sv=svd(X);

for i=1:L
    s(i)=sqrt(sum(sv(i+1:end).^2));
end

val=0;
for j=[s(1), s(3), s(10)]
    val=val+1;
    Ya(:,val)=mcw(Y,L,froweights(L,size(X,2)),j);
end

plot(Y)
hold on
P=plot(Ya);
NameArray = {'LineStyle'};
ValueArray = {'--',':','-.'}';
set(P,NameArray,ValueArray);
hold off


