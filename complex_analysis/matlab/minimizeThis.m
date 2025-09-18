function y=minimizeThis(x)
% This is the functional that is minimized by the demo script trySolve.m
%
global RAWDATA

N=length(x)/2;
p=x(1:N)+i*x(N+1:2*N);
F=RAWDATA.fisher;
X=RAWDATA.x;

y=real(p'*F*p)-2*real(p'*X);

