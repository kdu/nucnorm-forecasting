function [norms] = tau_norm_equiv(Y, N, L, ranks, w)
%TAU_NORM_EQUIV Computes equivalent norms for the approximations
  X=hmat(Y,L);
  
  norms = zeros(length(ranks),1);
  for j=1:length(ranks)
    Xr = lra(X, ranks(j));
    yr = hankvec_avg(Xr);

    norms(j) = norm(sqrt(w(:)).*(yr(:)-Y(:)));
  end  
end

