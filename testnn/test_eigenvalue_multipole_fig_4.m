clear all
clc
addpath ..;


% number of observations
n = 13;
d = n-1;

max_r = 4; % floor((n)/2);

% number of random experiments
%M = 100;

M = 20;
h11 = 1;
h1s = 0.7:0.05:1.2;

% for saving the errors
threshold = 1e-4;

%frob_real = zeros(length(h11), max_r);
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
      
      prob_complex(r,k,m) = prob_complex(r,k,m) + double(norm(p_complex(pr.tts) - ph_complex(pr.tts), 'fro') < threshold); 
   %   [frob_real(i,r)] = max(frob_real(i,r), ...
   %   [frob_complex(i,r)] = max(frob_complex(i,r), ...
   %                          norm(p_complex(pr.tts) - ph_complex(pr.tts), 'fro'));
    end
    prob_complex(r,:,:) / M
  end
end  
end
  
% load test

prob_complex / M


prob_r1 = reshape(prob_complex(1,:,:)/M, [length(h1s), d]);
f1 = figure('rend','painters','pos',[10 10 256 192])
%image(prob_r1 .* 255)

hh = pcolor(padarray(prob_r1, [1 1], 1, 'post'));
shading flat;
colormap('gray');
caxis([0 1]);
%grid off;
%set(hh, 'EdgeColor', 'none');
set(gca, 'YDir', 'normal')
%colormap('gray');
hold on

set(gca,'xTick',(1:d)+0.5)
set(gca,'xTickLabels',(1:d))
set(gca,'yTick',(1:2:length(h1s))+0.5)
set(gca,'yTickLabels', {'0.7', '0.8', '0.9', '1.0', '1.1', '1.2'})
y = 1-h1s.^(-2*(d+1))
mkr = log((sqrt(y.^2+4) + y) / 2) ./ log(abs(h1s))
mkr(isnan(mkr)) = n;
mkr = ceil(mkr) - 2;

xs = zeros(1,2*length(h1s));
ys = zeros(1,2*length(h1s));

ys(1:2:end) = (1:length(h1s));%-0.5;
ys(2:2:end) = (1:length(h1s))+1;%+0.5;
xs(1:2:end) = mkr+1;%0.5;
xs(2:2:end) = mkr+1;%0.5;

plot(xs, ys, 'r-', 'Linewidth', 3)

xlabel('m')
ylabel('rho')
export_fig_eps_own('m_vs_rho_r1.eps');

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

