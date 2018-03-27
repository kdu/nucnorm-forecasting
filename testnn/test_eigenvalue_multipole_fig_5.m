clear all
clc
addpath ..;


% number of observations
n = 13;
d = n-1;

max_r = floor((n)/2);

M = 20;
h11 = 1;
h1s = [0.85,0.95,1,1.05,1.15];

% for saving the errors
threshold = 1e-4;

prob_complex = zeros(max_r,length(h1s),d);

pr = buildSLRMCHankel(zeros(2*d+1, 1), d);


for r = 1:max_r

for k=1:length(h1s)
  k
  for j = 1:M
    j
    h11 = h1s(k);
    h_real = h11 * ones(1,r); %[(2* rand(1,r-1) - 1) 1]
    h_complex = h_real .* exp(1i * rand(1, r) * 2 * pi())
    
 
    p_complex = sum_1d_exp(h_complex, 2*n -1);
    pr_complex = pr; pr_complex.p(1:2*d+1) = p_complex(1:2*d+1);
    for m=1:d
      pr_complex.p(2*d+2-m) = NaN;
      
      ph_complex = nnSLRMC(pr_complex);
      
      prob_complex(r,k,m) = prob_complex(r,k,m) + double(norm(p_complex(pr.tts) - ph_complex(pr.tts), 'fro') < threshold); 
    end
    prob_complex(r,:,:) / M
  end
end  
end
  
load test5 % results were saved to test5.mat

for j=1:length(h1s)
  prob_h1 = reshape(prob_complex(:,j,:)/M, [max_r, d]);
  f1 = figure('rend','painters','pos',[10 10 256 192])

   hh = pcolor(padarray(prob_h1, [1 1], 1, 'post'));
   
   shading flat;
    colormap('gray');
    caxis([0 1]);
%  grid off;
%  set(hh, 'EdgeColor', 'none');
  set(gca, 'YDir', 'normal')


  set(gca,'xTick',(1:d)+0.5)
  set(gca,'xTickLabel',1:d)
  set(gca,'yTick',(1:max_r)+0.5)
  set(gca,'yTickLabel',1:max_r)

  
  
  xlabel('m')
  ylabel('r')
  
   export_fig_eps_own(sprintf('m_vs_r_h%1d_%02d.eps', floor(h1s(j)), mod(round(h1s(j)*100),100) ));
end
%  save2pdf('m_vs_rho_r1.pdf', f1);

% plot
% fontsize = 16;

% f1 = figure;
% colormap([0 0 0; 1 0 0])
% mesh(1:max_r,h11,frob_real,double(frob_real> 1e-5), 'LineWidth', 1.5);
% xlabel('rank','fontsize', fontsize);
% ylabel('rho','fontsize', fontsize);
% zlabel('Fro-error','fontsize', fontsize);
% title('Real roots, n = 9','fontsize', fontsize);
% save2pdf('hankel_n9_real.pdf', f1);
% 
% f2 = figure;
% colormap([0 0 0; 1 0 0])
% mesh(1:max_r,h11,frob_complex,double(frob_complex> 1e-5), 'LineWidth', 1.5);
% xlabel('rank','fontsize', fontsize);
% ylabel('rho','fontsize', fontsize);
% zlabel('Fro-error','fontsize', fontsize);
% title('Complex roots, n = 9','fontsize', fontsize);
% save2pdf('hankel_n9_complex.pdf', f2);

