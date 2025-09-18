#!/bin/bash



### BASIC DIRECTORIES ###########
BASE_DIR=../
DATA_DIR=$BASE_DIR/MPTA_4.5yr-example/
PSR_DIR=$DATA_DIR/data/

OUT_DIR=$BASE_DIR/MPTA_4.5yr-example/out/
SLURM_DIR=$OUT_DIR/slurm_out/
#################################


for lmax in 16;
	do
		nside=16
		cut=32
		binn=1

		level=0.5
		res=1

		NAME=MPTA-CGW


		# with calculations
		sbatch -p short.q --time=01:30:00 --mem=20GB --job-name=$NAME --output $SLURM_DIR/$NAME\_level$level\_res$res\_$lmax\_$nside\_$cut\_$n\.out --error $SLURM_DIR/$NAME\_level$level\_res$res\_$lmax\_$nside\_$cut\_$n\.err --wrap="singularity exec -B /scratch/ /scratch/kgrunthal/MPTA_singularity/ python3 $BASE_DIR/sensitivity_investigation.py --head_dir $OUT_DIR --psrpickle $PSR_DIR/psrs_sim.pkl --psrlist $PSR_DIR/psr_names.txt --tres $PSR_DIR/tres_MPTA.json --lmax $lmax --nside $nside --cutoff $cut --fbin $n --res $res --level $level --incCrossCorr --calculate"


		# without calculations 
		#sbatch -p short.q --time=01:30:00 --mem=20GB --job-name=$NAME --output $SLURM_DIR/$NAME\_level$level\_res$res\_$lmax\_$nside\_$cut\_$n\.out --error $SLURM_DIR/$NAME\_level$level\_res$res\_$lmax\_$nside\_$cut\_$n\.err --wrap="singularity exec -B /scratch/ /scratch/kgrunthal/MPTA_singularity/ python3 $BASE_DIR/sensitivity_investigation.py --head_dir $OUT_DIR --psrpickle $PSR_DIR/psrs_sim.pkl --psrlist $PSR_DIR/psr_names.txt --tres $PSR_DIR/tres_MPTA.json --lmax $lmax --nside $nside --cutoff $cut --fbin $n --res $res --level $level --incCrossCorr"
	done