%Vectors
function Ya=mc(Y,L,e)
N=length(Y);

cvx_begin sdp;
variable Yapp(N);
minimize(norm_nuc(hankel(Yapp(1:L),Yapp(L:N))))
subject to
    norm(hankel(Yapp(1:L),Yapp(L:N))-hmat(Y,L),'fro')<=e;
cvx_end
Ya=Yapp;

