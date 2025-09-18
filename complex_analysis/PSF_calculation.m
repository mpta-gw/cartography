function PSF_calculation(nside_in, lmax_in, sv_in)
% function PSF_calculation
%
% load clean map and covariance matrix
% -> transform both into pixel maps
%
% Careful: different lm order in MatLab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath("<path_to_matlab_folder>");
basepath = "<path_to_output_from_anisotropy_analysis>";

% choose "cc" or "nocc" depending on the analysis setting
wch = "cc";


sv = string(sv_in);
lmax = string(lmax_in);
nside = string(nside_in);
    
% loop over the frequency bins
for f_bin=1
  fbin = string(f_bin);
  disp(f_bin);

  % load (inverse) Fisher matrix and Covariance matrix
  
  M_pta = readmatrix(basepath+"/M_sph/M_sph_" + wch + "_" + lmax + "_" + sv + "_" + nside + "_fbin" + fbin + ".dat", "OutputType","double");
  Mprime_inv_pta = readmatrix(basepath+"/Mprimeinv_sph/Mprimeinv_sph_" + wch + "_" + lmax + "_" + sv + "_" + nside + "_fbin" + fbin + ".dat", "OutputType","double");
  pCov_pta = readmatrix(basepath+"/Cov_sph/Cov_sph_" + wch + "_" + lmax + "_" + sv + "_" + nside + "_fbin" + fbin + ".dat", "OutputType","double");
  

  % variables for matlab reordering
  M = 0*M_pta;
  pCov = 0*pCov_pta;
  Mprime_inv = 0*Mprime_inv_pta;

  % rearrange entries in order to match LIGO convention
  % - N = number of SpH modes [N = (l_max+1)^2]
  % - Lmax = maximum degree of SpH
  % - (lvec, mvec) = (l,m) order assumed by LIGO code
  % - (lvec_pta, mvec_pta) = (l,m) order assumed by PTA code

  N = length(M_pta(:,1));
  
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
    for jj=1:N
      idx_jj = lvec_pta==lvec(jj) & mvec_pta==mvec(jj);
      M(ii,jj) = M_pta(idx_ii, idx_jj);
      Mprime_inv(ii,jj) = Mprime_inv_pta(idx_ii, idx_jj);
      pCov(ii,jj) = pCov_pta(idx_ii, idx_jj);
    end
  end
  


  % pixel map resolution in degrees
  res=1;

  [~, ~, ~, ~, ~, U] = getSigmaMap(pCov, res);
  
  %=====================================================
  % collect all maps for saving to file
  writematrix(real(U), basepath+'/conversion_matrix/U-re_lmax'+lmax+'_res'+string(res) + '.txt', Delimiter='\t')
  writematrix(imag(U), basepath+'/conversion_matrix/U-im_lmax'+lmax+'_res'+string(res) + '.txt', Delimiter='\t')
  writematrix(real(Mprime_inv*M), basepath+'/maps/MpinvM-re_'+wch+'_lmax'+lmax+ '_res' + string(res) + '_sv'+string(sv)+'_fbin'+fbin + '.txt', Delimiter='\t')
  writematrix(imag(Mprime_inv*M), basepath+'/maps/MpinvM-im_'+wch+'_lmax'+lmax+ '_res' + string(res) + '_sv'+string(sv)+'_fbin'+fbin + '.txt', Delimiter='\t')
  
end

return
