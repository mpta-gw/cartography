function [lmcosi] = plm2plmreal (plm)

%  Takes input plm as complex spherical harmonics packed as:

%              |  plm
%   ---------------------------------
%   M=-L,L=Lmax|
%   .......    |
%   M=-1,L=Lmax|
%   .......    |
%   M=-1,L=1   |
%   M=0,L=0    |
%   M=0,L=1    |
%   .......    |
%   M=0,L=Lmax |
%   M=1,L=1    |
%   .......    |
%   M=1,L=Lmax |
%   .......    |
%   M=L=Lmax   |
%
%   and converts them to lmcosi matrix listing l,m,cosine and sine
%   coefficients of real spherical harmonics, in order of l, then m
%
%  Routine written by Madeleine Udell.
%  Contact madeleine_udell@yale.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lmax=sqrt(length(plm))-1;
sym=(length(plm)-(lmax+1))/2;

l=zeros((lmax+1)*(lmax+2)/2,1);
m=zeros((lmax+1)*(lmax+2)/2,1);
co=zeros((lmax+1)*(lmax+2)/2,1);
si=zeros((lmax+1)*(lmax+2)/2,1);



for ii=0:lmax;
    l((ii)*(ii+1)/2 + 1:(ii+1)*(ii+2)/2)=ii;
    m((ii)*(ii+1)/2 + 1:(ii+1)*(ii+2)/2)=(0:ii);
    co((ii)*(ii+1)/2 + 1)=real(plm(sym+1+ii));
    si((ii)*(ii+1)/2 + 1)=0; % has to be zero (S.Ballmer) imag(plm(sym+1+ii));
    % here, we treat ii as an index through the ms and jj for the ls
    if ii~=0;
    for jj=0:lmax-ii;
        co((ii+1)*(ii+2)/2 + jj*(jj+1)/2 + ii*jj) = real(plm((lmax-ii+1)*(lmax-ii+2)/2 - jj) + ((-1)^ii) * plm(length(plm) - (lmax-ii+1)*(lmax-ii+2)/2 + jj + 1))/sqrt(2);
        si((ii+1)*(ii+2)/2 + jj*(jj+1)/2 + ii*jj) = imag(plm((lmax-ii+1)*(lmax-ii+2)/2 - jj) - ((-1)^ii) * plm(length(plm) - (lmax-ii+1)*(lmax-ii+2)/2 + jj + 1))/sqrt(2);
	% Sign of si changed by Stefan Ballmer, 11/14/2007. This makes the transformatin consistent with the complex Ylm definition in Jackson
	% and the real spherical harmonics definition in the spherelib and LIGO-T070045-00-U.
    end;end;
end;  

lmcosi=[l,m,co,si];
