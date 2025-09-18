function [Cl,sig] = constructCl(plm, invfisher)
% Eric Thrane

% construct an array of Cl sigmavalues defined by complex plms:
%   
% C_l = (1/(2*l + 1)) * [sum_m |plm|^2 - |InvFisher_mm|]
%
% var_Cl = (2/(2*l + 1)^2) * sum_mm' |InvFisher_mm'|^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lmax = sqrt(length(plm))-1;

[lvec, mvec]=getLMvec(Lmax);

for l=0:Lmax
  ind = find(lvec==l);
  Cl(l+1) = sum(abs(plm(ind)).^2) - trace(abs(invfisher(ind,ind)));
  var(l+1) = sum(sum(abs(invfisher(ind,ind)).^2));        % from Joe's elog

  Cl(l+1) = Cl(l+1) / (2*l+1);
  sig(l+1) = sqrt(var(l+1)) * sqrt(2/(2*l+1)^2);
end

return
