function y = subX(x, Lsub)
%
% find the subset of the cross-correlation X_lm (packed in the 
% usual way in terms of l,m values) having l<=Lsub
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract Lmax from input array
N = size(x,1);
Lmax = sqrt(N) -1;
if Lsub>Lmax
  error('requesting subset bigger than original array');
end

% construct l,m vectors appropriate for this Lmax
[lvec, mvec] = getLMvec(Lmax);

% construct subset of the X_lm array having l<=Lsub
ind = find(lvec<=Lsub);
y = x(ind);

return
