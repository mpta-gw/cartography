#!/bin/bash

BASE_DIR=../
OUT_DIR=$BASE_DIR/MPTA_4.5yr-example/out/
SLURM_DIR=$OUT_DIR/slurm_out/

MATLAB_DIR=$BASEDIR/complex_analysis/matlab/

module load matlab/R2023a

sbatch -p short.q --time=00:10:00 --mem=10GB --output=$SLURM_DIR/matlab.out --error=$SLURM_DIR/matlab.err --job-name=ml_$i -D $MATLAB_DIR --wrap="matlab -nojvm -nodisplay -nosplash -nodesktop -r \"PSF_calculation(16, 8, 32)\""


