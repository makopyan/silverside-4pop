#!/bin/bash
export OMP_NUM_THREADS=9

BAMLIST=$1 # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
POP=$2
REFERENCE=$3 #Path to input saf file
OUTDIR=$4 # Path to output files
BASENAME=$5 # Basename for output files
NIND=$6
MIND=$7
MAXD=$8

#estimate the saf for all sites in genome:
/workdir/programs/angsd0.931/angsd/angsd -P 9 -b $BAMLIST -anc $REFERENCE -out $OUTDIR$POP$BASENAME -minMapQ 20 -minQ 20 -minInd $NIND -setMinDepth $MIND -setMaxDepth $MAXD -doCounts 1 -GL 1 -doSaf 1 -fold 1

#estimate the sfs:
/workdir/programs/angsd0.931/angsd/misc/realSFS $OUTDIR$POP$BASENAME'.saf.idx' -P 9 > $OUTDIR$POP$BASENAME'.sfs'


#how to run:
#nohup sh ./angsd_sfs_allsites.sh /workdir/arne/sampleinfo/BamFiles_clean_list_JIGA.txt JIGA /workdir/arne/ReferenceSeq/Mmenidia_refgenome_anchored.all_renamed_v2.fasta /workdir/arne/results/summary_statistics/ _allsites_folded 48 16 200 > jiga_unfolded_sfs_allsites_nohup.log &
