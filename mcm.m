%Matrices
function Ya=mcm(Y,L,e)
X=hmat(Y,L)
cvx_begin sdp;
variable Xa(L,length(Y)-L+1) hankel;
minimize(norm_nuc(Xa))
subject to
    norm(Xa-X,'fro')<=e;
cvx_end
Ya=[Xa(1,1:size(Xa,2)),Xa(2:end,end)']';