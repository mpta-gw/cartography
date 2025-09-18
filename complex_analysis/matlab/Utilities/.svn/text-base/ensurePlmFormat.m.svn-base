function plm=ensurePlmFormat(plm)

% input:  spherical harmonics coefficients, either complex or real format
%
% output: spherical harmonics coefficients, complex format and
%         ensuring that is represents a real distribution,
%         by projecting on to the real part.
%         

s=size(plm);

if s(2)==1 && mod(sqrt(s(1)),1)==0
  fprintf('\nInput is vector of complex spherical harmonics coefficients.\n\n');
  
elseif (s(2)==2 || s(2)==4 ) && mod((sqrt(8*s(1)+1)-3)/2,1)==0
  fprintf('\nInput is vector of real spherical harmonics coefficients.\n\n');
  plm=plmreal2plm(plm);
else
  error('No idea what the input is...');
end
 
s=size(plm);
lmax=sqrt(s(1))-1;
[l,m]=getLMvec(lmax);
plm=(plm + conj(mTranspose(plm)))/2;
