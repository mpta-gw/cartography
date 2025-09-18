# cartography
Codebase for nanohertz GWB anisotropy analysis and calculating PTA point-spread-function sky maps

## Usage
This code works on a three-step basis, as it contains both `python` and `matlab` scripts. The `python` scripts require a series of inputs, given via flags. In its current state, the framework also needs a distinct folder setup of the output directory of the analysis. Both the usage of the flags, as well as directory set up are encoded in the `bash` scripts provided in 'sbatch_scripts`.

1. ##### Setting up the PTA model and calculating the anisotropy analysis data products
    Calculate the anisotropy analysis products such as $M$, $X$, $\tilde{M}^{-1}$, $\mathcal{P}'$, etc. (see [Grunthal & Nathan et al. 2025](https://academic.oup.com/mnras/article/536/2/1501/7912549)) with the `python` script `OS_Anisotropy_Analysis_defiant.py`. It uses the regularised spherical harmonics anisotropy analysis. 
    The implementation makes use of the package [`defiant`](https://github.com/GersbachKa/defiant) to calculate the frequency-resolved pulsar pair correlations. The anisotropy analysis uses a modified version of the scripts from the [`MAPS`](https://github.com/NihanPol/MAPS)  package , altered to use complex spherical harmonics and providing maps in the observer-based GW polarisation basis, and customised for regular used of the regularisation scheme. 

    The analysis can be perfomed with or without regularisation, and with or without the full cross-correlation covariance matrix.

    The usage of `OS_Anisotropy_Analysis_defiant.py` is demonstrated in `sbatch_scripts/submit_OS-aniso.sh`.

  
2. ##### Converting the products to the pixel basis and calculate sky maps
    Conversion to the pixel basis is performed using `matlab`-based code developed by the LIGO collaboration (https://git.ligo.org/stochastic-public/stochastic).
    The `matlab` framework is hosted in  `complex_analysis/matlab/`, the individual scripts can be found in `complex_analysis`. It is necessary to explicitly include the path of this folder in the matlab script, otherwise it will not be able to call any function from the code framework.

    * `map_conversion_nomonopole.m` calculates a file containing the clean, clean S/N, radiometer S/N and sensitivity maps, and two other files with the maximum and minimum clean and radiometer S/N distributions. All of them are calculated with a removed monopole, ($\mathcal{P}'_{00} = 0$).
    * `map_conversion_nocleandistribution.m` calculates a file containing the clean, clean S/N, radiometer S/N and sensitivity maps, and one other file with the maximum and minimum radiometer S/N distributions. 
    * `PSF_calculation.m` produces the components of the PTA distortion matrix $\Lambda_{s_\mathrm{reg}}$ (cf. Grunthal et al. 2025b). These are two files with the real and imaginary part of $M \cdot \tilde{M}^{-1}$, and another two files with the real and imaginary part of the conversion matrix $U$. 


3. ##### Calculating the point-spread-function and plotting the result
    The evaluation of the PSF is implemented in `PSF_calculation_plotting.py`. 
    It loads the previously calculated matrices, and calculates the PSF by iterating through the columns of $\Lambda_{s_\mathrm{reg}}$. The PSF map is plotted, and also stored as a file. This allows replotting without recalculating.

    The usage of `PSFanalysis_plotting.py` is demonstrated in `sbatch_scripts/submit_PSFanalysis-plotting.sh`.



## Citation

This code is released as part of [], so if you use it this code, please consider referencing (and reading :D) this paper . 

It also
 - makes direct use of [`defiant`](https://github.com/GersbachKa/defiant), based on https://ui.adsabs.harvard.edu/abs/2025PhRvD.111b3027G/abstract
 - makes direct use of [stochastic GW analysis code](https://git.ligo.org/stochastic-public/stochastic) developed by the LIGO collaboration
 - uses a modified version of [`MAPS`](https://github.com/NihanPol/MAPS)

