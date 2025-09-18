function retval = checkPixelConversion(Lmax,res)
% Verifies that the global variable
% PIXELCONVERSION is set for the appropriate Lmax and res values
%
% * If not: try to load the file sigmaMapData_r_l.mat, which is located
%   in the LEGENDRE directory of the spherelib package.
%   Here r is the resoluiton rounded to 1/100 of a degree and l = Lmax
%
% * If that file does not exist: calculate the conversion data and save the
%   file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%   Lmax           : Lmax desired for the following analysis
%   res            : Resolution res desired for the following analysis
%
% Output:
%   retval         : =1 if PIXELCONVERSION was set correctly
%                  | =2 if PIXELCONVERSION was set with loaded data
%                  | =3 if PIXELCONVERSION was set with calculated data
%                  | =0 if an error occured
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The global variable PIXELCONVERSION has the following fields:
%   Lmax           : Lmax - this is verified
%   res            : resolution res - this is verified
%   U              : conversion matrix that takes (complex) SpH basis
%                  | to pixel basis, i.e. X_pix = U * X_sph
%   RA             : Vector of right ascentions
%   DECL           : Vector of declinations
%   dOmg           : Vector of steradians covered by each pixel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global PIXELCONVERSION;

retval=0;

try
    checkFlag=and(PIXELCONVERSION.Lmax==Lmax,PIXELCONVERSION.res==res);
catch
    checkFlag=false;
end

if checkFlag
    retval=1;    
else
    % check whether the required file exists
    resstr=num2str(res*100,'%5.0f');
    lmaxstr=num2str(Lmax,'%5.0f');
    pathstr = fileparts(mfilename('fullpath'));
%EHT: on 240212
%EHT    LEGENDRE = [pathstr '/../lib/spherelib/LEGENDRE'];
    LEGENDRE = ['./spherelib/LEGENDRE'];
    file=fullfile(getenv('IFILES'),...
        [LEGENDRE,'/pixelConversionData_',resstr,'_',lmaxstr,'.mat']);
    if exist(file,'file')
        disp(['Loading conversion data for Lmax= ',...
               num2str(Lmax),', res= ',num2str(res),...
               ' from file ',file]);
        load(file);
        if and(PIXELCONVERSION.Lmax==Lmax,PIXELCONVERSION.res==res)
            retval=2;
        else
            retval=0;
            error(['Data file ',file,' not consistent']);
        end
    else
        disp(['Calculating conversion data for Lmax= ',...
               num2str(Lmax),', res= ',num2str(res)]);
        %Prepare the structure
        dum.Lmax=Lmax;
        dum.res=res;
        PIXELCONVERSION=dum;
        Nra=ceil(360/res);
        Ndecl=ceil(180/res)+1;
        N=(Lmax+1)^2;
        PIXELCONVERSION.U    = zeros(Nra*Ndecl,N);
        %calculate the data
        ylm=eye(N);
        [PIXELCONVERSION.U(1:Nra*Ndecl,1),...
            PIXELCONVERSION.RA,PIXELCONVERSION.DECL,...
            PIXELCONVERSION.dOmg] = makemap(ylm(:,1),res,0,1);
        for ii=2:N
            PIXELCONVERSION.U(1:Nra*Ndecl,ii) = makemap(ylm(:,ii),res,0,1);
        end
        retval=3;
        % save the file
        disp(['Saving conversion data for Lmax= ',...
               num2str(Lmax),', res= ',num2str(res),...
               ' to file ',file]);
        save(file,'PIXELCONVERSION');
    end
end

end

