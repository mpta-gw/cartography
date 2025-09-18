function [p, fisher, x, invfisher] = invertX(filename);

% Give the filename of the index file of the results from testSpH_____.m
index = load (filename);

% The number of iterations of the code
imax = length(index.Sky);

% Initialize variables
loadeddatafile='';
dummy=load (index.Sky{1}.filename);
x=zeros(size(dummy.Sky{1}.X));
fisher=zeros(size(dummy.Sky{1}.Fisher));

% Sum the projection X and the Fisher matrix for each iteration
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
end;

invfisher = inv(fisher);
p = invfisher * x;

