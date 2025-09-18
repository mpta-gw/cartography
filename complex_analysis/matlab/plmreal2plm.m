function [plm,lvec,mvec] = plmreal2plm(lmcosi)

% converts the real sperical harmonics format lmcosi to
% the complex sperical harmonics format plm used by the
% stochastic SpH analysis

s=size(lmcosi);

% calculate lmax from number of rows
lmax=(sqrt(8*s(1)+1)-3)/2;

if mod(lmax,1)~=0
 error('Number of rows must be consistent the a lmax (acceptable: 1,3,6,10, etc)');
 plm=-1;
 return;
end


if s(2)==4
  data = lmcosi(:,3:4);
  lm   = lmcosi(:,1:2);
else
  data = lmcosi;
  lm   = zeros(s(1),2);
  ii=1;
  for ll=0:lmax
    for mm=0:ll
      lm(ii,1:2)=[ll,mm];
      ii=ii+1;
    end
  end
end

plm=zeros((lmax+1)^2,1);
lvec=plm;
mvec=plm;

for ii=1:s(1)
  l=lm(ii,1);
  m=lm(ii,2);
  jp=getIndex( m,l,lmax);
  if m==0
    plm(jp,1) =transpose(data(ii,1));
    lvec(jp,1)=l;
    mvec(jp,1)=m;
  else
    jm=getIndex(-m,l,lmax);
    msgn=(-1).^m;
    % I wounld have chosen this matrix (consistent with LIGO-T070045-00-U.)
    M=[msgn,-i*msgn;1,i]/sqrt(2);
  % Madeleine's choice
  % M=[msgn,i*msgn;1,-i]/sqrt(2);
    plm([jp,jm],1) =M*transpose(data(ii,:));
    lvec([jp,jm],1)=l;
    mvec([jp,jm],1)=[m;-m];
  end
end


% test code
% mmm=[0];
% lll=[0];
% ind=getIndex(mmm,lll,0);
% [mmm,lll,ind]
% mmm=[-1;0;0;1];
% lll=[1;0;1;1];
% ind=getIndex(mmm,lll,1);
% [mmm,lll,ind]
% mmm=[-2;-1;-1;0;0;0;1;1;2];
% lll=[2;2;1;0;1;2;1;2;2];
% ind=getIndex(mmm,lll,2);
% [mmm,lll,ind]
% mmm=[-3;-2;-2;-1;-1;-1;0;0;0;0;1;1;1;2;2;3];
% lll=[3;3;2;3;2;1;0;1;2;3;1;2;3;2;3;3];
% ind=getIndex(mmm,lll,3);
% [mmm,lll,ind]


return

function ind=getIndex(m,l,lmax)
  ind=1 + lmax*(lmax+1)/2 + m.*(2*lmax+1-abs(m))/2 + l.*(m>=0) + (lmax-l).*(m<0);  
return
