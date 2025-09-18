function [dems,dels,mzero,lmcosi,mzin,mzout]=addmon(L)
% [dems,dels,mzero,lmcosi,mzin,mzout]=ADDMON(L)
%
% Returns a vector with spherical harmonic orders m up until
% a given degree 'L'. For real spherical harmonic orders, m=0:L
% The idea is to make a matrix LMCOSI as input/output for
% programs XYZ2PLM and PLM2XYZ.
%
% INPUT:
%
% L           Degree of the expansion
%
% OUTPUT:
%
% dems        Lists all orders m involved in the expansion
% dels        Lists all degrees l involved in the expansion
% mzero       Index to the m=0 locations for the zonal coefficients
% lmcosi      Matrix ready to put coefficients in
% mzin        Location of m=0 sin terms in vector not listing them
% mzout       Running index into lmcosi(3:4) to return vector without sin(0)
%
% Last modified by fjsimons-at-alum.mit.edu, May 3rd, 2004

worist=[addmup(0:L,'a')+1];
iench=[0 ; ones(worist(end)-1,1)];
wadist=-(0:L);
iench(worist)=wadist;
dems=cumsum(iench(1:end-1));
dels=gamini(0:L,[0:L]+1);
dels=dels(:);
mzero=[1 ; worist(1:end-1)'];
lmcosi=[dels dems zeros(length(dels),2)];
mzin=(mzero*2+1)-[1:length(mzero)]';

indo=reshape(1:length(dems)*2,length(dems),2);
indo(length(dems)+mzero)=0; indo=indo'; indo=indo(:);
mzout=indo(indo~=0);

% diff(sort(mzout)) should be 1, 2 or 3
