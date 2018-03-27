M=6;
N=72;
L=24;

load deathsdata.mat
Yfull=deathsdata;
Y=deathsdata(1:N);

trs=[];
for tau=1000:1000:15000
for rho=0.01:0.01:0.25
    for i=1:N
        w(i)=exp(rho*i);
    end
Yapp=mcwf(Y,L,M,w',tau);


smse=sqrt((1/M)*sum((Yapp(N+1:N+M)-Yfull(N+1:N+M)).^2));
trs=[trs; tau, rho, smse];
end
end

surfc(1000:1000:15000, 0.01:0.01:0.25, log10(reshape(trs(:,3), 25, 15)))
xlabel('tau')
ylabel('alpha')
zlabel('log10(rootMSE)')
min(trs(:,3))
