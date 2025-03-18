#!/bin/bash
export OMP_NUM_THREADS=9

BAMLIST=$1 # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
POP=$2
REFERENCE=$3 #Path to input saf file
OUTDIR=$4 # Path to output files
BASENAME=$5 # Basename for output files
REFINDEX=$6
WINDOWS=$7

#estimate the saf for all sites in genome:
/workdir/programs/angsd0.931/angsd/angsd -P 18 -b $BAMLIST -anc $REFERENCE -fai $REFINDEX -rf $WINDOWS -out $OUTDIR$POP$BASENAME -minMapQ 20 -minQ 20 -GL 1 -doSaf 1
