function [p_res] = nnSLRMC(slra_pr, is_quiet)
%
% For complex case: the direct usage of cvx with complex variables gives
% poor performance. 
% We have to construct an equivalent real structure instead (its description
% can be found in the appendix of the paper).
%
  if nargin < 2
    is_quiet = true;
  end    
  tts = slra_pr.tts;
  p0 = slra_pr.p;
  p0 = p0(:);
  nm_ind =  find(~isnan(p0));
  
  
  
  if isreal(p0(nm_ind))
    cvx_begin, cvx_quiet(is_quiet);
        variable ph(length(p0))
      minimize norm_nuc(ph(tts))
      subject to
        ph(nm_ind) == p0(nm_ind);
    cvx_end;
  else
    len = length(p0);
    tts2 = [tts       (tts+2*len); ...
            (tts+len) (tts+3*len)];
    
    p0r = real(p0(nm_ind));
    p0i = imag(p0(nm_ind));
    cvx_begin, cvx_quiet(is_quiet);
        variable ph4(4*len)
      minimize norm_nuc(ph4(tts2))
      subject to
        ph4(nm_ind) == p0r;
        ph4(len+nm_ind) == p0i;
        ph4(1:len) == ph4((1:len) + 3*len);
        ph4((1:len)+len) == -ph4((1:len) + 2*len);
    cvx_end;    
    
    ph = ph4(1:len) + 1i* ph4((1:len)+len);
  end    

  
  p_res = full(ph);
end