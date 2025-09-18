function x=sbank(c,h0,f0,Pa,Ua,Kp,Ku,lev,n,an,dn,pph,intel)
% x=SBANK(c,h0,f0,P,U,Kp,Ku,lev,n,an,dn)
% x=SBANK(c,h0,f0,P,U,Kp,Ku,lev,n,an,dn,pph,intel)
%
% Filterbank: SYNTHESIS, one iteration, one branch
% See: ABANK, FBANK, WT, IWT
%
% Note: this function is invoked iteratively in IWT;
% We're calculating projections on V 
% (for the scaling coefficient) and on W
% (for the wavelet coefficients);
% We don't add them yet to find the multiresolution
% projections; the successive addition is done by 
% IWT or SUCCAPP.
%
% 'pph' 1= Time-domain full bitrate (inefficient);
%       2= Time-domain polyphase (inefficient);  
%       3= Z-domain polyphase (default)
%       4= Lifting
% 'intel' 1 Integer operations
%         0 Full precision
%
% 'c'        wavelet or scaling coefficients
% 'h0','f0'  the filter coefficients
% 'P','U','Kp','Ku' the lifting operators
% 'lev'      the last branch they've been on:
%            'a' approximation/lowpass
%            'd' detail/highpass
% 'n'        the number of cascades the coefficients have experienced
% 'an','dn'  the number of coefficients in each branh at each level
%            so we can figure out how many coefficients to reconstruct.
%
% This is a recursive algorithm.
% See example at IWT.
% Assume the lengths of the arrays after the entire reconstruction is EVEN

% Last modified by fjsimons-at-alum.mit.edu, 3/18/2003

% Only max(an) is used to figure out what the last level of recursion is

defval('pph',3)
defval('intel',0)
defval('deb',0) % Debugging output

% Make conjugate mirror filters from the factors of the product filter
[h0,h1,f0,f1,l,lf0,lh0,a0,a1]=prodco(f0,h0);

% Synthesis bank
switch lev
  % Which filter generated this level?
 case 'a' 
  lzf=length(h0);
  dv= 'd2hx';
  % Work with F0 all the way back down
  zfil=f0;
  al= 'a0';
  alp= +1;
  if pph==4
    y0=c/Kp;
    y1=repmat(0,dn(n),1)/Ku;
  end
 case 'd'
  lzf=length(f0);
  dv= 'd2fx';
  % Work with F1 for this step n
  zfil=f1;
  al= 'a1';
  alp= -1;
  if pph==4
    y0=repmat(0,an(n),1)/Kp;
    y1=c/Ku;
  end
end

% Number of coefficients of the signal you're reconstructing
% At the last level you need to assume previous level was EVEN
% Otherwise, previous level might have been odd, so take that one.
% But if the filter itself is of ODD length, need to add one again

if pph~=4
  if n>1
    lx=max(an(n-1),2*length(c)-lzf);
  else
    lx=2*length(c)-lzf+mod(lzf,2);
  end
  [d2fx,d2hx]=landd(f0,h0,1:lx);
  dv=eval(dv);
else
  if n>1
    lx=an(n-1);
  else
    lx=2*length(c);
  end
end

% Get downsampling arrays; use for upsampling
switch pph
 case 4
  if n==1 & lev=='d' ; disp('Lifting implementation'); end
  % UNDO SCALE ------------------------------------
  M=1;
  if iscell(Pa); M=length(Pa); end

  % Reverse loop over lifting steps-------------------
   for index=M:-1:1
     if iscell(Pa)
       P=Pa{index};
       U=Ua{index};
     else
       P=Pa;
       U=Ua;
     end
     Pl=length(P);
     Ul=length(U);

     % UNDO UPDATE -----------------------------------
     for l=1+floor(Ul/2):dn(n)-ceil(Ul/2)+1
       Lu=l-[floor(Ul/2):-1:(1-ceil(Ul/2))];
       if intel==1; y0(l)=y0(l)-floor(U(:)'*y1(Lu)+1/2); end
       if intel==0; y0(l)=y0(l)-U(:)'*y1(Lu); end
     end
     if deb==1
       disp(sprintf(['Undo Update  %i /%i with ',...
		     repmat('%i ',size(Lu)),'/%i'],...
		    l,length(y0),Lu,length(y1)))
     end
     
     % UNDO PREDICT ----------------------------------
     % The last condition only when the prediction operator
     % is of length unity and when the a's and d's have
     % unequal length
     for l=ceil(Pl/2):an(n)-floor(Pl/2)-(Pl==1)*(an(n)~=dn(n))
       Lp=l+[(1-ceil(Pl/2)):floor(Pl/2)];
       if intel==1; y1(l)=y1(l)+floor(P(:)'*y0(Lp)+1/2); end
       if intel==0; y1(l)=y1(l)+P(:)'*y0(Lp); end
     end  
     if deb==1
       disp(sprintf(['Undo Predict %i /%i with ',...
		     repmat('%i ',size(Lp)),'/%i'],...
		    l,length(y1),Lp,length(y0)))
     end
   end
     
   % MERGE -----------------------------------------
   x(~even(1:lx),:)=y0;
   x(even(1:lx),:)=y1;

 case 3
  if n==1 & lev=='d' ; disp('Z-domain polyphase implementation'); end
  % Make this half Polyphase - we're doing branch by branch
  % Here's what we'll do differently next: full polyphase
  al=eval(al);
  ZFeven=zfil(~~(al+alp));
  ZFodd=zfil(~(al+alp));
  lcd=logical(1/2-1/2*(-1).^(1:length(dv)+length(zfil)-1))';
  lzfec=length(ZFeven)+length(c);
  lzfoc=length(ZFodd)+length(c);

  y0=repmat(0,floor(length(lcd)/2),1);
  y1=repmat(0,ceil(length(lcd)/2),1);
  % xout=repmat(NaN,length(lcd),1);
  y0(1:lzfec-1)=conv(ZFeven(:),c(:));
  y1(1:lzfoc-1)=conv(ZFodd(:),c(:));

  % Even output
  x(~lcd,1)=y0;
  % Odd output
  x(lcd,1)=y1;  
  % Undo delay corresponding to the odd degree of the product filter
  x=x(l+1:end-l);
 case 2
  if n==1 & lev=='d' ; disp('Time-domain polyphase implementation'); end
  % Synthesis
  % Size of F0 is (lx+lh0-1+lf0-1) by (lx+lh0-1)
  % Size of F1 is (lx+lf0-1+lh0-1) by (lx+lf0-1)
  % and element (1,1) is h0(1)
  % Filter Toeplitz matrices
  ZF=convmtx(zfil(:),length(dv));
  % Construct the polyphase SYNTHESIS matrix
  % Polyphase form
  % Upsample: delete columns: keep these
  ZF=ZF(:,dv);
  % This already reduces the number of multiplications by two
  % Now split in odd and even
  xout=repmat(NaN,length(dv)+length(zfil)-1,1);  
  ZFeven=ZF(even(xout),:);
  ZFodd=ZF(~even(xout),:);
  PS=[ZFeven ; ZFodd];
  % Could make this TYPE II (p. 132)
  % The output of this is mixed and it needs to be recombined!  
  % Check out the operation PS*PA*X - you've switched evens 
  % and odds and added a delay here and there so the 
  % reconstruction automatically makes up for it
  xoutm=PS*c;
  % Output of the original evens of x
  xout(even(xout),1)=xoutm(1:size(ZFeven,1));
  % Output of the original odds of x
  xout(~even(xout),1)=xoutm(size(ZFeven,1)+1:end);
  x=xout;
  % Undo delay corresponding to the odd degree of the product filter
  x=x(l+1:end-l);
 case 1
  if n==1 & lev=='d' ; disp('Toeplitz implementation'); end
  % Or else be inefficient
  % Upsample
  blk=~(1:length(dv));
  blk(dv)=c;
  % Filter
  ZF=convmtx(zfil(:),length(dv));
  x=ZF*blk(:);
    % Undo delay corresponding to the odd degree of the product filter
  x=x(l+1:end-l);
end

% Recursive algorithm - always go back on the lowpass tree
if n>1
  x=sbank(x,h0,f0,Pa,Ua,Kp,Ku,'a',n-1,an,dn,pph,intel);
end
