clear all
clc
addpath ..;


% number of observations
n = 13;
d = n-1;

max_r = floor((n)/2);

% number of random experiments
%M = 100;

M = 20;
h11 = 1;

% for saving the errors
threshold = 1e-4;

%frob_real = zeros(length(h11), max_r);
prob_complex = zeros(max_r,d);

pr = buildSLRMCHankel(zeros(2*d+1, 1), d);

for r=1:max_r
  r
  for j = 1:M
    j
    h_real = h11 * ones(1,r); %[(2* rand(1,r-1) - 1) 1]
    h_complex = h_real .* exp(1i * rand(1, r) * 2 * pi())
    
%     figure;
%     plot(complex(exp(1i * linspace(0,2*pi(),200))));
%     hold on
%     plot(h_complex, '*');
    
    p_complex = sum_1d_exp(h_complex, 2*n -1);
    pr_complex = pr; pr_complex.p(1:2*d+1) = p_complex(1:2*d+1);
    for m=1:d
%      m
      pr_complex.p(2*d+2-m) = NaN;
      
      ph_complex = nnSLRMC(pr_complex);
      
      prob_complex(r,m) = prob_complex(r,m) + double(norm(p_complex(pr.tts) - ph_complex(pr.tts), 'fro') < threshold); 
   %   [frob_real(i,r)] = max(frob_real(i,r), ...
   %   [frob_complex(i,r)] = max(frob_complex(i,r), ...
   %                          norm(p_complex(pr.tts) - ph_complex(pr.tts), 'fro'));
    end
    prob_complex / M
  end
end  
  

prob_complex / M
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

