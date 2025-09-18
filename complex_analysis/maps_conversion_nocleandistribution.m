function maps_conversion_nocleandistribution(basepath_in, nside_in, lmax_in, sv_in)
% function maps_coversion_nocleandistribution
%
% load clean map and covariance matrix
% -> transform both into pixel maps
%
% Careful: different lm order in MatLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath("<path_to_matlab_folder>");
basepath = string(basepath_in);

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

  % clean map, inverse Fisher matrix and Covariance matrix for S/N map
  pOpt_pta = readmatrix(basepath+"/popt_sph/popt_sph_cc_"+lmax+"_"+sv+"_"+nside+"_fbin"+fbin+".dat", "OutputType","double");
  pOpt_pta = pOpt_pta(:,2);
  pCov_pta = readmatrix(basepath+"/Cov_sph/Cov_sph_cc_"+lmax+"_"+sv+"_"+nside+"_fbin"+fbin+".dat", "OutputType","double");
  Mprime_inv_pta = readmatrix(basepath+"/Mprimeinv_sph/Mprimeinv_sph_cc_"+lmax+"_"+sv+"_"+nside+"_fbin"+fbin+".dat", "OutputType","double");

  %save('mpta_complex.mat', 'pOpt_pta', 'pCov_pta');

  % variables for matlab reordering
  X = 0*X_pta;
  M = 0*M_pta;
  pOpt = 0*pOpt_pta;
  pCov = 0*pCov_pta;
  Mprime_inv = 0*Mprime_inv_pta;

  % rearrange entries in order to match LIGO convention
  % - N = number of SpH modes [N = (l_max+1)^2]
  % - Lmax = maximum degree of SpH
  % - (lvec, mvec) = (l,m) order assumed by LIGO code
  % - (lvec_pta, mvec_pta) = (l,m) order assumed by PTA code

  N = length(pOpt);
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
    pOpt(ii) = pOpt_pta(idx_ii);
    for jj=1:N
      idx_jj = lvec_pta==lvec(jj) & mvec_pta==mvec(jj);
      M(ii,jj) = M_pta(idx_ii, idx_jj);
      Mprime_inv(ii,jj) = Mprime_inv_pta(idx_ii, idx_jj);
      pCov(ii,jj) = pCov_pta(idx_ii, idx_jj);
    end
  end
  


  % pixel map resolution in degrees
  res=1;

  %-- 1. Radiometer map --------------------------------
  [diag_fisher, ~, ~, ~, U] = diagPixel(M, res);
  X_pix = U*X;
  p_radio = real(X_pix./diag_fisher);
  s_radio = real(diag_fisher.^-0.5);
  
  sn_radiometer = p_radio./s_radio;
 

  %-- 2.  map -------------------------------------
  [clean_map, ~, ~] = makemap(pOpt, res);
  

  %-- 3. S/N map ---------------------------------------
  % - sigmaPix:     the sigma map
  % - (RA, DEC):    arrays of RA and DEC for mapping

  [~, sigmaPix, RA, DEC, ~, U] = getSigmaMap(pCov, res);
  sn = clean_map./sigmaPix;

  %-- 4. sensitivity map -------------------------------
  Lambda = U*Mprime_inv*M*U';
  sn_sensitivity = real(diag(Lambda))./sigmaPix;


  %=====================================================
  % collect all maps for saving to file
  writematrix([RA, DEC, clean_map, sigmaPix, sn, sn_radiometer, sn_sensitivity], basepath+'/maps/maps_'+lmax+'_'+sv+'_'+nside+'_fbin'+fbin+'.txt', Delimiter='\t');
  writematrix([RA, DEC, real(diag_fisher)], basepath+'/maps/map_M_'+lmax+'_'+sv+'_'+nside+'_fbin'+fbin+'.txt', Delimiter='\t');
 
  
  % estimate significance of maximum and minimum in radiometer map
  
  [V, D] = eig(M);
  ev_M = diag(D);
  Trials = 1e4;
  %f = waitbar(0,'Trials: 0/10000');

  for kk=1:Trials
    %text='Trials: %d/10000';
    %waitbar(kk/1e4, f, sprintf(text, kk));
    
    % create random vector of X 
    %  - eigenvalues are distributed with sigma=sqrt(ev_M)
    %  - calculate X using V from the random evs
    ev_Msqrt_rand = sqrt(real(ev_M)).*randn(N, 1);
    X_rand = V*ev_Msqrt_rand;

    % construct radiometer S/N map from drawn X
    %  - convert X using transformation matrix U
    %  - divide by fisher matrix in pix basis
    % dirty map in pixel basis
    X_pix_rand = U*X_rand;
    p_radio_rand = real(X_pix_rand./diag_fisher);


    % record max and min snr
    radio_max_snr(kk) = max(p_radio_rand./s_radio);
    radio_min_snr(kk) = min(p_radio_rand./s_radio);
  end


  % collect max and min vectors for saving to file
  writematrix([transpose(radio_max_snr),transpose(radio_min_snr)], basepath+'/distribution/distribution_radiometer_maxmin_'+lmax+'_'+sv+'_'+nside+'_fbin'+fbin+'.txt', Delimiter='\t');

  

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
