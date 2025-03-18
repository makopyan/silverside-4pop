#!/bin/bash

BAMLIST=$1 # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
POP=$2
REFERENCE=$3 #Path to input saf file
OUTDIR=$4 # Path to output files
BASENAME=$5 # Basename for output files
NIND=$6
MIND=$7
MAXD=$8
WINDOWS=$9

#estimate the saf for all sites in genome:
/workdir/programs/angsd0.931/angsd/angsd -P 18 -b $BAMLIST -anc $REFERENCE -rf $WINDOWS -out $OUTDIR$POP$BASENAME -minMapQ 20 -minQ 20 -minInd $NIND -setMinDepth $MIND -setMaxDepth $MAXD -doCounts 1 -GL 1 -doSaf 1 -fold 1
