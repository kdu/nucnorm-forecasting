clear
%set(0,'defaulttextinterpreter','latex')
yf=[];

L = 15;
%L = 25;
K = 50;

rho = 0.01:0.01:0.2;

for n = 80:20:140
    y=[];
for x = rho
    y = [y,Cfun(x,n,L,K)];
end
    yf=[yf;y];
end
figure
plot(rho,yf)
h=legend('$n=80$', '$n=100$','$n=120$','$n=140$')

set(h,'Interpreter','latex');
legend('Location','best')
legend('boxoff')