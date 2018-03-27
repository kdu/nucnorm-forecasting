%Weighted vectors
function Ya=mcw(Y,L,w,e)
N=length(Y);

cvx_begin sdp;
%cvx_solver mosek;
%cvx_precision low;
variable Yapp(N);
minimize(norm_nuc(hankel(Yapp(1:L),Yapp(L:N))))
subject to
    norm(sqrt(w).*(Y-Yapp))<=e;
cvx_end
Ya=Yapp;

