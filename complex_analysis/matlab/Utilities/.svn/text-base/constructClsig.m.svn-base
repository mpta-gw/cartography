function sigCl = constructClsig(plm, invfisher)
% construct an array of Cl sigmavalues defined by complex plms:
%   
% C_l = (1/(2*l + 1)) * sum_m=-l^l abs(plm)^2
%
% sigma_Cl = (4/(2*l + 1)) * sqrt( sum abs(plm)^2 sigma_lm^2 )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lmax = sqrt(length(plm))-1;

lmvec = 1:1:(1+Lmax)^2;
[lvec, mvec]=getLMvec(Lmax);

sigCl = zeros(Lmax+1,1);

for l=0:Lmax
  ind = find(lvec==l);
  for m=1:length(ind)
    for n=1:length(ind)
%Small signal approximation--------------------------------------
%      sigCl(l+1) = sigCl(l+1) + abs(invfisher(ind(m),ind(n)))^2;
%----------------------------------------------------------------
      sigCl(l+1) = sigCl(l+1) + ...
      abs(invfisher(ind(m),ind(n)))^2 + ...
      abs( conj(plm(ind(m)))*invfisher(ind(m),ind(n))*plm(ind(n)) ) + ...
      abs(conj( conj(plm(ind(m)))*invfisher(ind(m),ind(n))*plm(ind(n)) ));
%WRONG-----------------------------------------------------------
%    plm(ind(m))*invfisher(ind(m),ind(n))*plm(ind(n))  + ...
%    ctranspose( (plm(ind(m))*invfisher(ind(m),ind(n))*plm(ind(n)) ));
    end
  end
  L=l-1;
  sigCl(l+1) = sqrt(sigCl(l+1)) * ( sqrt(2)/(2*L+1) );
end

return
