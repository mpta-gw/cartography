function addo=adrc(mat,varargin)
% addo=ADRC(mat)
% addo=ADRC(mat,value)
%
% Adds row and column to 'mat' for use wih PCOLOR
%
% Last modified by fjsimons-at-alum.mit.edu, October 15th, 2003

if nargin==1
  val=0;
else
  val=varargin{1};
end

[m,n]=size(mat);
addo=[mat  repmat(val,m,1) ; repmat(val,1,n+1)];
