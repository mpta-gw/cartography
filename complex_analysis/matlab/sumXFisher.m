function [p, x, fisher, invfisher, coh, chisquared, Psi] = sumXFisher(filenames);

for kk=1:length(filenames)
filename=filenames{kk};
% Give the filename of the index file of the results from testSpH_____.m
index = load (filename);

% The number of iterations of the code
imax = length(index.Sky);

% Initialize variables
existPsi = 1;
loadeddatafile='';
dummy=load (index.Sky{1}.filename);
if kk==1
  x=zeros(size(dummy.Sky{1}.X));
  fisher=zeros(size(dummy.Sky{1}.Fisher));
  coh=zeros(size(dummy.Sky{1}.coh));
  try 
    Psi=zeros(size(dummy.Sky{1}.Psi));
  catch
    Psi=0;
    existPsi=0;
  end;
end

% Sum the projection X, the Fisher matrix for each iteration
for ii=1:imax
    datafilename = index.Sky{ii}.filename;
    if strcmp(datafilename,loadeddatafile);
        jj=jj+1; % an index to keep track of where we are in a file
    else
        loadeddatafile=datafilename;
        datafile = load (datafilename);
        jj=1;
    end;
    x = x + datafile.Sky{jj}.X;
    fisher = fisher + datafile.Sky{jj}.Fisher;
    coh = coh + datafile.Sky{jj}.coh;
    if existPsi
      Psi = Psi + datafile.Sky{jj}.Psi;
    end
end;

end;

invfisher = reginv(fisher);
p = invfisher * x;

% take real part to avoid imaginary contributions from ill-conditioning
chisquared = coh - ctranspose(p)*x - ctranspose(x)*p + real(ctranspose(p)*fisher*p);

return

