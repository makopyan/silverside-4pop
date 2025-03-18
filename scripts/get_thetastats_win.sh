#!/bin/bash
export OMP_NUM_THREADS=9

BAMLIST=$1 # Path to textfile listing bamfiles to include in global SNP calling with absolute paths
POP=$2
REFERENCE=$3 # Path to reference genome
SAFFILE=$4 #Path to input saf file
OUTDIR=$5 # Path to output files
BASENAME=$6 # Basename for output files
WIN=$7 #Window size in bp for sliding window analysis
STEP=$8 #Step size in bp for sliding window analysis, set equal to WIN for non-overlapping windows

#estimate the sfs used as prior
/workdir/programs/angsd0.931/angsd/misc/realSFS $SAFFILE -P 4 > $OUTDIR$POP$BASENAME'.sfs'

#estimate per-site theta
/workdir/programs/angsd0.931/angsd/angsd -bam $BAMLIST -out $OUTDIR$POP$BASENAME -doThetas 1 -doSaf 1 -pest $OUTDIR$POP$BASENAME'.sfs' -anc $REFERENCE -GL 2 -fold 1

#estimate genome-wide neutrality stats
/workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat $OUTDIR$POP$BASENAME'.thetas.idx'

#estimate window-based neutrality stats
/workdir/programs/angsd0.931/angsd/misc/thetaStat do_stat $OUTDIR$POP$BASENAME'.thetas.idx' -win $WIN -step $STEP -outnames $OUTDIR$POP$BASENAME'_theta.thetasWindow.gz'


#Example for running script:
#nohup sh ./angsd_thetastats_win.sh /workdir/arne/sampleinfo/BamFiles_clean_list_JIGA.txt JIGA /workdir/arne/ReferenceSeq/Mmenidia_refgenome_anchored.all_renamed_v2.fasta /workdir/arne/results/summary_statistics/JIGA_filtsnps.saf.idx /workdir/arne/results/summary_statistics/ _thetastats 10000 10000 > /workdir/arne/output_logfiles/jiga_thetastats_folded_nohup.log &  
