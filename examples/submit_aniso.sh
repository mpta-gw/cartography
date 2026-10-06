#!/bin/bash
date +"%T"

PTA=pta_1

### RELEVANT DIRECTORIES ############################################

# head path
HP=/u/kgrunthal/IPTA/ipta_aniso_datachallenge/

# working directory
WD=$HP/scripts/

# data directory
DD=$HP/datasets/$PTA

# output directory
OUTDIR=$HP/outdir/$PTA

#####################################################################


### ENSURE EXISTENCE OF OUTPUT FOLDERS ##############################
if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR/
fi


if [ ! -d "$OUTDIR/figures/" ]; then
    mkdir $OUTDIR/figures/
fi

#####################################################################


### SETTINGS ANISOTROPY ANALYSIS ####################################
lmax=29
nside=16
cut=40
#####################################################################



for fbin in 1;
    do
        NAME=aniso_$lmax\_$cut\_$fbin

        # no pickle file available
        #singularity exec -B /u/ /u/kgrunthal/environments/MPTA_singularity/ python3 $WD/OS_Anisotropy_Analysis_with_conversion.py --datapath $DD --psrlist $HP/datasets/psrlist_full.txt --noisefile $WD/dictionaries/WNnom_GWB.json --noisemodels $WD/dictionaries/WN_CRN.json    --fbin $fbin --lmax $lmax --nside $nside --cutoff $cut --incCrossCorr      --basepath $WD --create_plots --outdir $OUTDIR

        # pulsar pickle available
        singularity exec -B /u/ /u/kgrunthal/environments/MPTA_singularity/ python3 $WD/OS_Anisotropy_Analysis_with_conversion.py --datapath $DD --psrlist $HP/datasets/psrlist_full.txt --psrpickle $DD/psrs.pkl --noisefile $WD/dictionaries/WNnom_GWBnoA.json --noisemodels $WD/dictionaries/WN_CRN.json    --fbin $fbin --lmax $lmax --nside $nside --cutoff $cut --incCrossCorr      --basepath $WD --create_plots --outdir $OUTDIR


    done

date +"%T"
