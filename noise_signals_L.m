

clear

N=100;
b=1;

sigma=0.1;
N_sim=50;

w=ones(N-1,1);

Lind=0;
for L=[10,20,30,40,50];
Lind=Lind+1;
    for sim=1:N_sim;
    for n=1:N
        s(n)=(b^n)*cos(2*pi*n/10);
    end
    
    
    %r=sigma*randn(N,1);
    r=sigma*ar1( N, 0.5);
    
    
    Y=s+r';
    X=hmat(Y(1:end-1),L);
        
    sv=svd(X);
    for i=1:L
        s1(i)=sqrt(sum(sv(i+1:end).^2));
    end

    
    for k=2:2
        un(k)=sqrt(wvnorm(hankvec(lra(X,k))-Y(1:end-1)',w));
    end
    
    
    for i=2:2
    Ya(:,sim,i)=mcwf(Y(1:end-1)',L,1,froweights(L,size(X,2)),s1(i));
    Ya_unit(:,sim,i)=mcwf(Y(1:end-1)',L,1,w,un(i));

    rmse(sim,Lind,i)=sqrt((Ya(end,sim,i)-Y(end))^2);
    rmse_unit(sim,Lind,i)=sqrt((Ya_unit(end,sim,i)-Y(end))^2);
    end
end
end

   

for i=2:2
    figure
    boxplot(rmse(:,:,i))
    ax = gca;
    ax.XTickLabel = {'10','20','30','40','50'};
    ylim([0 0.4])
    %export_fig(sprintf('cosw%d.eps', i));
end

for i=2:2
    figure
    boxplot(rmse_unit(:,:,i))
    ax = gca;
    ax.XTickLabel = {'10','20','30','40','50'};
    ylim([0 0.4])
    %export_fig(sprintf('cosw_unit%d.eps', i));
end
    