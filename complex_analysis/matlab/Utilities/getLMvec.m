function [lvec,mvec]=getLMvec(Lmax)

for M=0:Lmax


  % index into the gammaLM matrix for the M=L values (M=0,1, ..., Lmax)
  i0 = M*(Lmax+1)- M*(M-1)/2 + 1;

  % vector of m values, l values
  mvec(i0:(i0+Lmax-M),1) = M;
  lvec(i0:(i0+Lmax-M),1) = M:Lmax;

end

% extend vectors of l and m values to negative values of m
lvec = [ flipud(lvec(2+Lmax:end)); lvec];
mvec = [-flipud(mvec(2+Lmax:end)); mvec];

