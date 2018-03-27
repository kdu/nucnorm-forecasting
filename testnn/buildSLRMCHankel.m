function s = buildSLRMCqHankel(p, d) 
  s.p = NaN * ones(1, 2*d+1);
  s.p(1:(d+1)) = p(1:(d+1));
  s.tts = hankel(1:(d+1), (d+1):(2*d+1));
end