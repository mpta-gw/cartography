function R=zRotation(Lmax,alpha,diagonly)

% Returns a rotation matrix for a rotation along the y axis
% for the complex spherical harmonics representation.
% This is the trivial rotation

% Since R is a diagonal Matrix, there is an option to just return the diagonal.

try, diagonly; catch diagonly=false; end;

[lvec,mvec]=getLMvec(Lmax);

% Sign: rotate coefficients in opposite direction
r=exp(-i*alpha*mvec);

if diagonly
  R=r;
else
  R=diag(r);
end;
