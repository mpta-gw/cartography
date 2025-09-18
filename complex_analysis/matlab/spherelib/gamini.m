function bigger=gamini(data,folding)
% bigger=GAMINI(data,folding)
% Replicates the ith element of a VECTOR 'folding'(i) times.
% If 'folding' is a scalar takes same for all elements.
%
% Example:
%
% [a,b]=degamini(gamini([1 2 3 1 4 5 ],[1 2 3 2 4 2]))
%
%
% See also GAMINI2
%
% Last modified by fjsimons-at-alum.mit.edu, June 8th 2000

if prod(size(folding))==1
  folding=repmat(folding,prod(size(data)),1);
end

data=data(:)';
folding=folding(:)';

if ~all(size(data)==size(folding)) 
  error([ ' Sizes of input and folding must be the same'])
end

gelp=zeros(1,sum(folding));
gelp([1 cumsum(folding(1:end-1))+1])=1; 
bigger=data(cumsum(gelp));
