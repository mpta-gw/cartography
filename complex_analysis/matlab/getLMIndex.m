function ind=getLMIndex(m,l,lmax)
  if l>lmax
    error('l > lmax');
  end
  if abs(m)>l
    error('abs(m) > l');
  end
  ind=1 + lmax*(lmax+1)/2 + m.*(2*lmax+1-abs(m))/2 + l.*(m>=0) + (lmax-l).*(m<0);
return
