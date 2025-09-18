% hack to enforce symmetry of Fisher matrix
fprintf('hack to enforce Fisher matrix symmetry!\n');
for ii=1:N
  for jj=1:ii-1
    fisher(ii,jj) = conj(fisher(jj,ii));
  end
end
	
% make a complex-version
if 1==0
fprintf('bug here?! hack to convert to complex sph\n');
X_tmp = X;
fisher_tmp = fisher;
for ii=1:N
  idx_ii = lvec==lvec(ii) & mvec==-mvec(ii);
  if mvec(ii)<0
    X(ii) = (1/sqrt(2))*(X_tmp(idx_ii) - i*X_tmp(ii));
	elseif mvec(ii)==0
	  X(ii) = X_tmp(ii);
  else
    X(ii) = ((-1)^ mvec(ii)/sqrt(2))*(X_tmp(idx_ii) + i*X_tmp(ii));
  end
  % update Fisher too?
	for jj=1:N
    idx_jj = lvec==lvec(jj) & mvec==-mvec(jj);
  end
end
end
	
