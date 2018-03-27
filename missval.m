
clear

N=100;


sigma=0.1;
N_sim=50;
L=30;
k=2;

    for n=1:N
        s(n)=cos(2*pi*n/10); % case 1
        %s(n)=cos(2*pi*n/10)*exp(0.02*n); % case 2
    end

for m=1:15;
    for sim=1:N_sim;
        r=sigma*randn(N,1);
        Y=s+r';
        X=hmat(Y(1:end-m),L);

        wF=froweights(L,size(X,2));
        wUnit=ones(N-m,1);
        
        clear wExp;
        for i=1:N-m
        wExp(i)=1.03^i;
        end
        wExp=wExp';
        
        tauF=sqrt(wvnorm(hankvec_avg(lra(X,k))-Y(1:end-m)',wF))
        tauUnit=sqrt(wvnorm(hankvec_avg(lra(X,k))-Y(1:end-m)',wUnit))
        tauExp=sqrt(wvnorm(hankvec_avg(lra(X,k))-Y(1:end-m)',wExp))
        
    Ya(:,sim,m)=mcwf(Y(1:end-m)',L,m,wF,tauF);
    Ya_unit(:,sim,m)=mcwf(Y(1:end-m)',L,m,wUnit,tauUnit);
    Ya_exp(:,sim,m)=mcwf(Y(1:end-m)',L,m,wExp,tauExp);

    rmse_F(sim,m)=sqrt((1/m)*sum((Ya(N-m:end,sim,m)-Y(N-m:end)').^2));
    rmse_unit(sim,m)=sqrt((1/m)*sum((Ya_unit(N-m:end,sim,m)-Y(N-m:end)').^2));
    rmse_exp(sim,m)=sqrt((1/m)*sum((Ya_exp(N-m:end,sim,m)-Y(N-m:end)').^2));
    end

end

%load missval_res1

figure('rend','painters','pos',[10 10 280 225])
boxplot(rmse_F)
axis([0.5 15.5 0 0.4]);
export_fig_eps_own(sprintf('rmse_F_1.eps', i));


figure('rend','painters','pos',[10 10 280 225])
boxplot(rmse_unit)
axis([0.5 15.5 0 0.4]);
export_fig_eps_own(sprintf('rmse_unit_1.eps', i));


figure('rend','painters','pos',[10 10 280 225])
boxplot(rmse_exp)
axis([0.5 15.5 0 0.4]);
export_fig_eps_own(sprintf('rmse_exp_1.eps', i)); 


% %load missval_res2
% 
% figure('rend','painters','pos',[10 10 280 225])
% boxplot(rmse_F)
% axis([0.5 15.5 0 0.6]);
% export_fig_eps_own(sprintf('rmse_F_2.eps', i));
% 
% 
% figure('rend','painters','pos',[10 10 280 225])
% boxplot(rmse_unit)
% axis([0.5 15.5 0 0.6]);
% export_fig_eps_own(sprintf('rmse_unit_2.eps', i));
% 
% 
% figure('rend','painters','pos',[10 10 280 225])
% boxplot(rmse_exp)
% axis([0.5 15.5 0 0.6]);
% export_fig_eps_own(sprintf('rmse_exp_2.eps', i)); 
    