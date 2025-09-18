function out = mTranspose(in, doRows, doColumns)

% calculates the (-1)^m*v_l-m for vectors and matrices

s=size(in);
lmax1=sqrt(s(1))-1;
lmax2=sqrt(s(2))-1;

try, doRows;    catch
  doRows   = not(mod(lmax1,1));
end;
try, doColumns; catch
  doColumns= not(mod(lmax2,1));
end;

if mod(lmax1,1)~=0 && doRows
  warning('Number of Rows not consistent with an lmax. No transposing along 1st dim.')
  doRows=false;
end

if mod(lmax2,1)~=0 && doColumns
  warning('Number of Columns not consistent with an lmax. No transposing along 2nd dim.')
  doColumns=false;
end

if s(1)>1 && s(2)>1
  if doColumns
    for ii=1:s(1)
      in(ii,:)=mTranspose(in(ii,:));
    end
  end
  if doRows
    for jj=1:s(2)
      in(:,jj)=mTranspose(in(:,jj));
    end
  end
  out=in;
else
  if s(2)>s(1)
    in=transpose(in);
  end
  % vector
  lmax=sqrt(length(in))-1;
  if mod(lmax,1)~=0
    error('Vector length not consistent with an lmax')
  end
  j0=1 + lmax*(lmax+1)/2;
  out=[flipud(in(j0+lmax+1:end));in(j0:j0+lmax);flipud(in(1:j0-1))];  
  [lv,mv]=getLMvec(lmax);
  out=out.*(-1).^mv;
  if s(2)>s(1)
    out=transpose(out);
  end
end


return
