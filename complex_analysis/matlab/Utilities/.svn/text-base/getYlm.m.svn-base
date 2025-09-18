function plm=getYlm(Lmax,l,m,cs);

% This function returns a vector corresponding to
% a basis vector in the complex spherical harmonics representation.
% input:
% Lmax: defines size of vector
% l:    l quantum number
% m:    m quantum number (either 0..l for real or -l..l for complex basis vector
% cs:   if omitted: complex Ylm,
%       if 0: cosine basis vector for l,m
%       if 1: sine  basis vector for l,m
%
% Written by Stefan Ballmer
%

if abs(m)>l
  error('Need -l <= m <= l');
end
if l>Lmax || abs(m)>Lmax
  error('Need l <= Lmax');
end

if exist('cs','var')
  if cs~=0 && cs~=1
    error('Third index cs has to be either 0 or 1.');
  end
  lmcosi=zeros((Lmax+1)*(Lmax+2)/2,4);
  ii=1;
  for ll=0:Lmax
    for mm=0:ll
      lmcosi(ii,1:2)=[ll,mm];
      if l==ll && m==mm
        lmcosi(ii,cs+3)=1;
      end
      ii=ii+1;
    end
  end
  plm=plmreal2plm(lmcosi);
else
  plm=zeros((Lmax+1)^2,1);
  plm(getIndex(m,l,Lmax),1)=1;
end;



return

function ind=getIndex(m,l,lmax)
  ind=1 + lmax*(lmax+1)/2 + m.*(2*lmax+1-abs(m))/2 + l.*(m>=0) + (lmax-l).*(m<0);  
return
