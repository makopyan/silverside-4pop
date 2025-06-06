#!/bin/bash
## This script is used to get minor allele frequency estimation from angsd for each population / group

BASEDIR=$1 #  Path to the base directory where adapter clipped fastq file are stored in a subdirectory titled "adapter_clipped" and into which output files will be written to separate subdirectories. An example for the Greenland cod data is: /workdir/cod/greenland-cod/
SAMPLETABLE=$2 # Path to a sample table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID, the 2nd column is the lane number, and the 3rd column is sequence ID. The combination of these three columns have to be unique. The 6th column should be data type, which is either pe or se. An example of such a sample table is: /workdir/cod/greenland-cod/sample_lists/sample_table.tsv
POPCOLUMN=$3 # The column index of the variable that you want to group by in the sample table above. In the Greenland project, it's the fifth column,  and thus 5
BAMLISTPREFIX=$4 # Prefix of the bam lists. An example from the Greenland cod project is bam_list_realigned_
REFERENCE=$5 # Path to reference genome
SNPLIST=$6 # Path to the SNP list
MINDP=$7 # Minimum depth filter
MAXDP=$8 # Maximum depth filter
MININD=$9 # Minimum individual filter
MINQ=${10} # Minimum quality filter
MINMAPQ=${11:-20} # Minimum mapping quality (alignment score) filter, default value is 20
ANGSD=${12:-/workdir/programs/angsd0.931/angsd/angsd} # Path to ANGSD, default value is /workdir/programs/angsd0.931/angsd/angsd
THREADS=${13:-28} # Number of parallel threads to use, default value is 28.
EXTRA_ARG=${14:-'-remove_bads 1 -only_proper_pairs 1 -C 50'} # Extra arguments when running ANGSD, default value is '-remove_bads 1 -only_proper_pairs 1 -C 50'

OUTBASE=`echo $SNPLIST | sed 's/\..*//' | sed -e 's#.*\/\(\)#\1#'`
CHRLIST=`echo $SNPLIST | sed 's/\..*//'`.chrs
OUTDIR=$BASEDIR'angsd/popminind'$MININD'/'
if [ ! -d "$OUTDIR" ]; then
	mkdir $OUTDIR
fi

for POP in `tail -n +2 $SAMPLETABLE | cut -f $POPCOLUMN | sort | uniq`; do 
	echo $POP
	$ANGSD \
	-b $BASEDIR'sample_lists/bam_list_per_pop/'$BAMLISTPREFIX$POP'.txt' \
	-anc $REFERENCE \
	-ref $REFERENCE \
	-out $OUTDIR$POP'_'$OUTBASE'_popminind'$MININD \
	-dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
	-P $THREADS \
	-setMinDepth $MINDP -setMaxDepth $MAXDP -minInd $MININD -minQ $MINQ -minMapQ $MINMAPQ \
	-sites $SNPLIST -rf $CHRLIST \
	$EXTRA_ARG \
	>& $BASEDIR'nohups/'$POP'_'$OUTBASE'_popminind'$MININD'_maf.log'
done
