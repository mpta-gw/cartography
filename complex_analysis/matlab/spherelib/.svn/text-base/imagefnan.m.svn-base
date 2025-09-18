function varargout=imagefnan(c11,cmn,matrix,colmap,cax,nancol,invo)
% IMAGEFNAN(c11,cmn,matrix,colmap,cax,nancol,invo)
% h=IMAGEFNAN
% 
% IMAGEF for NaN-filler and multiple colormaps
%
% See IMAGEF, IMAGEFDIR, IMAGENAN, JOINCOLMAP
%
% Last modified by fjsimons-at-alum.mit.edu, Jan 23th 2002

defval('nancol',[1 1 1])
defval('colmap',[kelicol])
defval('invo',0)
defval('c11',[0 1])
defval('cmn',[1 0])
defval('cax',[0 1])

if isstr(colmap)
  colormap(colmap)
  colmap=colormap;
  if invo==1
    colmap=flipud(colmap);
  end
end

ncol=linspace(cax(1),cax(2),size(colmap,1))';
thenans=find(isnan(matrix));

matrix(thenans)=0;
matrix(matrix<cax(1))=cax(1);
matrix(matrix>cax(2))=cax(2);

X=interp1(ncol,colmap,matrix(:),'nearest');

if ~isempty(thenans),
  X(thenans(:),1)=nancol(1);
  X(thenans(:),2)=nancol(2);
  X(thenans(:),3)=nancol(3);
end

X=reshape(X,[size(matrix),3]);

cm1=[c11(1) cmn(2)];
c1n=[cmn(1) c11(2)];

h=image([cm1(1) c1n(1)],[cm1(2) c1n(2)],flipdim(X,1)); 

axis xy image

if nargout
  varargout{1}=h;
end


  
