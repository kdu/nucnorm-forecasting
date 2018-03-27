%Weighted vectors with forecast
function Ya=mcwf(Y,L,M,w,e)
N=length(Y);

cvx_begin sdp;
variable Yapp(N+M);
minimize(norm_nuc(hankel(Yapp(1:L),Yapp(L:N+M))))
subject to
    norm(sqrt(w).*(Y(1:N)-Yapp(1:N)))<=e;
cvx_end
Ya=Yapp;
