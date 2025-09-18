#!/bin/bash


# folders
BASE_DIR=../
DATA_DIR=$BASE_DIR/MPTA_4.5yr-example/

OUT_DIR=$BASE_DIR/MPTA_4.5yr-example/out/
SLURM_DIR=$OUT_DIR/slurm_out/


### CREATING FOLDERS ############

if [ ! -d "$OUT_DIR" ]; then
    mkdir $OUT_DIR/
fi

if [ ! -d "$SLURM_DIR" ]; then
    mkdir $SLURM_DIR/
fi


if [ ! -d "$OUT_DIR/figures/" ]; then
    mkdir $OUT_DIR/figures/
    mkdir $OUT_DIR/figures/OS_Anisotropy_diagnostics/
fi

if [ ! -d "$OUT_DIR/Mprimeinv_sph/" ]; then
    mkdir $OUT_DIR/M_sph/
    mkdir $OUT_DIR/Mprimeinv_sph/
    mkdir $OUT_DIR/Cov_sph/
    mkdir $OUT_DIR/X_sph/
    mkdir $OUT_DIR/distribution/
    mkdir $OUT_DIR/maps/
    mkdir $OUT_DIR/popt_sph/
    mkdir $OUT_DIR/conversion_matrix
fi



# Analysis parameters
lmax=8
nside=16
cut=32


for bin in 1;
	do
        NAME=OSdef_AN--MPTA_CGW--$lmax\_$cut\_$n

       	sbatch -p short.q --time=02:00:00 --mem=0 --output=$SLURM_DIR/$NAME\.out --error=$SLURM_DIR/$NAME\.err --job-name=$NAME --wrap="singularity exec -B /scratch/ /scratch/kgrunthal/MPTA_singularity/ python3 $BASE_DIR/OS_Anisotropy_Analysis_defiant.py --data $DATA_DIR/data/ --psrlist $DATA_DIR/data/psr_names.txt --utilfolder $BASE_DIR/PTA_AnalysisUtils/ --SPNAmodels $DATA_DIR/noisefiles/MKT_SPNAmodels_WN.json -noisefile $DATA_DIR/noisefiles/MKT_noisevalues_simWNGWB.json --outdir $OUT_DIR/ --ptamodel $DATA_DIR/ptamodel-OS_CRN.json --fbin $bin --lmax $lmax --nside $nside --cutoff $cut --create_plots - --incCrossCorr"


    done




