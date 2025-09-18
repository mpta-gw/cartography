function R=WignerRotation(Lmax,beta)

% Returns a rotation matrix for a rotation along the y axis
% for the complex spherical harmonics representation
% This is the non-trivial rotation

[lvec,mvec]=getLMvec(Lmax);

T=(Lmax+1)^2;
R=zeros(T,T);

for l=0:Lmax
  ii=find(lvec==l);
  R(ii,ii)=WignerDMatrix(l,beta);
end

