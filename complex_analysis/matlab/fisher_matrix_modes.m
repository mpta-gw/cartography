function fisher_matrix_modes(basepath, nside_in, lmax_in, sv_in, mode)
% function fisher_matrix_modes
%
% load clean map and covariance matrix
% -> transform both into pixel maps
%
% Careful: different lm order in MatLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("/scratch/kgrunthal/");

sv = string(sv_in);
lmax = string(lmax_in);
nside = string(nside_in);

% loop over the first three frequency bins
for f_bin=1:3
  fbin = string(f_bin);
  disp(f_bin);
  % X and M for radiometer map
  X_pta = readmatrix(basepath+"/X_sph/X_sph_cc_"+lmax+"_"+sv+"_"+nside+"_fbin"+fbin+".dat", "OutputType","double");
  X_pta = X_pta(:,2);
  M_pta = readmatrix(basepath+"/M_sph/M_sph_cc_"+lmax+"_"+sv+"_"+nside+"_fbin"+fbin+".dat", "OutputType","double");


  %save('mpta_complex.mat', 'pOpt_pta', 'pCov_pta');

  % variables for matlab reordering
  X = 0*X_pta;
  M = 0*M_pta;

  % rearrange entries in order to match LIGO convention
  % - N = number of SpH modes [N = (l_max+1)^2]
  % - Lmax = maximum degree of SpH
  % - (lvec, mvec) = (l,m) order assumed by LIGO code
  % - (lvec_pta, mvec_pta) = (l,m) order assumed by PTA code

  N = length(X);
  Lmax = sqrt(N)-1;

  [lvec, mvec] = getLMvec(Lmax);

  lvec_pta = [];
  mvec_pta = [];
  for ll=0:Lmax
    for mm=-ll:ll
      lvec_pta = [lvec_pta ll];
      mvec_pta = [mvec_pta mm];
    end
  end

  % carry out conversion
  for ii=1:N
    idx_ii = lvec_pta==lvec(ii) & mvec_pta==mvec(ii);
    X(ii) = X_pta(idx_ii);
    for jj=1:N
      idx_jj = lvec_pta==lvec(jj) & mvec_pta==mvec(jj);
      M(ii,jj) = M_pta(idx_ii, idx_jj);
    end
  end


  % pixel map resolution in degrees
  res=1;

  % calculate radiometer map
  [U,S,V] = svd(M);


  Pn = 0*X;
  Pn(mode) = 1;
  popt = U*Pn;

  [cleanmap, RA, DEC] = makemap(popt, res);

  
  writematrix([RA, DEC, cleanmap], basepath+'/modes/mode'+string(mode)+'_lmax'+lmax+'_nside'+nside+'_fbin'+fbin+'.txt', Delimiter='\t');

end

return

% COMMENT THE FOLLOWING IN
% IF YOU WANT TO PLOT WITH MATLAB

% make plots of everything... starting with SNR
%  the -1 number is for a variable called spt, which has a default value = 0.
%  I think you can use this pass clean and sigma separately.
plotMapAitoff(clean_map./sigmaPix, 360, 181, -1);
print('-dpng', '/scratch/kgrunthal/MeerKAT_PTA/complex_analysis/img/snr')

% plot sigma map
plotMapAitoff(sigmaPix, 360, 181, -1);
print('-dpng', '/scratch/kgrunthal/MeerKAT_PTA/complex_analysis/img/sigma')

% plot clean map
plotMapAitoff(clean_map, 360, 181, -1);
print('-dpng', '/scratch/kgrunthal/MeerKAT_PTA/complex_analysis/img/clean')

return
