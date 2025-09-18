function val=fralmanac(neem,plenet)
% val=FRALMANAC('name','planet')
%
% Accesses a database with planetary constants given in SI units.
% If 'name' is empty returns a list of possibilities.
%
% INPUT: 
%
% 'name'     'DegDis'    length of one equatorial degree of longitude [m]
%            'Radius'    Volumetric mean radius                       [m]
%            'GravCst'   Gravitational constant                       [m^3kg^{-1}s^{-2}]
%            'GravAcc'   Surface gravity                              [ms^{-2}]
%            'IMR2'      Reduced moment of inertia (I/MR2)            [dimensionless]
%            'a_EGM96'   Earth Reference radius                       [m]
%            'GM_EGM96'  Earth Reference mass constant                [m^3s^{-2}]
%            'a_GLGM2'   Mars Reference radius                        [m]
%            'GM_GLGM2'  Mars Reference mass constant                 [m^3s^{-2}]
%            'a_GLGM2'   Moon Reference radius                        [m]
%            'GM_GLGM2'  Moon Reference mass constant                 [m^3s^{-2}]
%
%            'EGM96'     Earth Geopotential 
%            'GMM2B'     Mars Geopotential (from MGS)
%            'JGM85H02'  Mars Geopotential (from MGS)
%            'GLTM2B'    Moon Topography (from Clementine)
%            'GTM3AR'    Earth Topography (from Georg Wenzel)
%            'Mars2000'  Mars Topography (from MOLA)
%            'GTM090'    Mars Topography (from MOLA)
%
% 'planet'   'Earth' (default)
%            'Mars'
%            'Moon'
%            'SHM' spherical harmonics models
%            'XYZ' spatially expanded models
%
% Amongt others, from http://nssdc.gsfc.nasa.gov/planetary/factsheet/
%
% If you are just copying this file only, you must supply the data in
% your own directories... 
%
% Last modified by fjsimons-at-alum.mit.edu, Jan 21th, 2004

defval('neem',[])
defval('plenet','Earth')

load(fullfile(getenv('IFILES'),'EARTHMODELS','CONSTANTS',plenet))

if nargin>0 & ~isempty(neem) 
    val=eval([plenet '.' neem]);
    if iscell(val) & prod(size(val))==1
      val=val{1};
    end
else
  number(str2mat(fieldnames(eval(plenet))))
end




