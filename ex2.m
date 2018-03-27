%ex2 deaths data with forecasting
M=6;
N=72;
L=24;

load deathsdata.mat
Yfull=deathsdata;
Y=deathsdata(1:N);

X=hmat(Y,L);
Xsize=hmat(Yfull,L);
w1=froweights(L,size(X,2));

sv=svd(X);

for i=1:L
    s(i)=sqrt(sum(sv(i+1:end).^2));
end

w2 = ones(1,N);
w3 = zeros(1,N);
lambda=0.05;
for i=1:N
    w3(i)=exp(lambda*i);
end

ranks = [3,6,12];

tau_w1 = s(ranks);
% tau_w1 = tau_norm_equiv(Y, N, L, ranks, w1); 
% equal to s(ranks) if hankvec_avg is used
tau_w2 = tau_norm_equiv(Y, N, L, ranks, w2);
tau_w3 = tau_norm_equiv(Y, N, L, ranks, w3);

s(ranks)
tau_w1
tau_w2
tau_w3

Ya = zeros(length(ranks)*3,N+M); 

for i=1:length(ranks), 
  Ya(i+0*length(ranks),:) = mcwf(Y,L,M,w1(:),tau_w1(i));
  Ya(i+1*length(ranks),:) = mcwf(Y,L,M,w2(:),tau_w2(i)); 
  Ya(i+2*length(ranks),:) = mcwf(Y,L,M,w3(:),tau_w3(i)); 
end

smse = zeros(size(Ya,1),1);


for i=1:size(Ya,1)
  figure('rend','painters','pos',[10 10 256 192])
  plot(Yfull);
  hold on
  plot(Ya(i,:),'k--');
  line([72 72], [min(Y(:))-10 max(Y(:))+10 ]);
  axis([0 length(Yfull) min(Y(:))-20 max(Y(:))+20] );
  hold off
  
  export_fig_eps_own(sprintf('f%dnew.eps', i));
  smse(i)=sqrt((1/M)*sum((Ya(i,N+1:N+M)'-Yfull(N+1:N+M)).^2));
end  


% Create a cell array of strings
j = 1;
str = {};
for i=1:length(ranks)
  str = [str;'$W_1$, rank = ',num2str(ranks(i)),sprintf(' & %6.0f',Ya(j,N+1:end)),...
      sprintf('& %6.2f\\\\', smse(j))];
  j = j+1;
end
for i=1:length(ranks)
  str = [str;'$W_2$, rank = ',num2str(ranks(i)),sprintf(' & %6.0f',Ya(j,N+1:end)),...
      sprintf('& %6.2f\\\\', smse(j))];
  j = j+1;
end
for i=1:length(ranks)
  str = [str;'$W_3$, rank = ',num2str(ranks(i)),sprintf(' & %6.0f',Ya(j,N+1:end)),...
      sprintf('& %6.2f\\\\', smse(j))];
  j = j+1;
end
% Save to file
fid = fopen('nn_forecasts.txt', 'w');
fprintf(fid, '%s\n', str{:});
fclose(fid);



