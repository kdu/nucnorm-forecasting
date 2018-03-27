N = length(fort);

x1 = hankel(Ya(1:L,1),Ya(L:N,1));
x2 = hankel(Ya(1:L,2),Ya(L:N,2));
x3 = hankel(Ya(1:L,3),Ya(L:N,3));

s1=svd(x1)
s2=svd(x2)
s3=svd(x3)

% figure
% hist(sqrt(s1(1:15)))
% figure
% hist(sqrt(s2(1:15)))
% figure
% hist(sqrt(s3(1:15)))

figure
plot(sqrt(s1(1:15)))
figure
plot(sqrt(s2(1:15)))
figure
plot(sqrt(s3(1:15)))