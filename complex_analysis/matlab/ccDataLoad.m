function [list,x]=ccDataLoad(file)
% function [list,x]=ccDataLoad(file)
%
%  loads a ccstats time series or ccstats frequency series data
%  file and converts it into a cell array.
%  
%  input:  file: filename 
%
%  output: list: cell array with data
%          x:    vector of time shifts or frequencies
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  table=load(file);
catch
  fprintf('Unable to load file %s\n',file);
  list=[];
  t=[];
  return;
end;
if length(table)==0
  fprintf('File %s is empty\n',file);
  list=[];
  t=[];
  return;
end;
l=size(table,1);

  l=length(table(:,1));
  Nsegs=length(unique(table(:,1)));
  N=l/Nsegs;

cols=size(table,2);
switch cols,
  case 3,
    x=0; % time shift zero
  case 4,
    x=table(1:N,3);
  case 5,
    x=table(1:N,3);
    if or(any(x<0),any(x>8192))
      error('Unknown file format: inconsistend x axis data');
    end
  otherwise,
    error('Unknown file format');
end;
  
for kk=0:(Nsegs-1)
  list{kk+1}.time=table(kk*N+1,1);
  switch cols,
    case 3,
      list{kk+1}.data=table(kk*N+1:kk*N+N,[2,3]); 
      list{kk+1}.ccstat=1;
    case 4,
      % real data
      list{kk+1}.data=table(kk*N+1:kk*N+N,[4,2]); 
      list{kk+1}.ccstatsTimeSeries=1;
    case 5,
      % complex data
      list{kk+1}.data(1:N,1)=table(kk*N+1:kk*N+N,4)+i*table(kk*N+1:kk*N+N,5); 
      list{kk+1}.data(1:N,2)=table(kk*N+1:kk*N+N,2); 
      list{kk+1}.ccspec=1;
    otherwise,
      error('Unknown file format');
    end;
end;
return;
