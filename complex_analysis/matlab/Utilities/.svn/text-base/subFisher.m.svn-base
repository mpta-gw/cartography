function y = subFisher(fisher, Lsub)
%
% find the subset of a Fisher information matrix (packed in the 
% usual way in terms of l,m values) having l<=Lsub
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract Lmax from input Fisher matrix
N = size(fisher,1);
Lmax = sqrt(N) -1;
if Lsub>Lmax
  error('requesting subset bigger than original fisher matrix');
end

% construct l,m vectors appropriate for this Lmax
[lvec, mvec] = getLMvec(Lmax);

% construct subset of fisher matrix having l<=Lsub
ind = find(lvec<=Lsub);
y = fisher(ind,ind);

return
