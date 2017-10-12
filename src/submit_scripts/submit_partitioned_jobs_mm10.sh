#! /bin/bash

# submit partitioned specificity score jobs
indir=$1

jobs=$(ls -1 $indir/Cas9*.txt)

for j in $jobs
do
	# make output directory
	[[ $j =~ ^.*\_([0-9]+).txt ]]
	partition=${BASH_REMATCH[1]}
        outdir="$(dirname $indir)/specificity_output_vs_mm10/partition_$partition"

	# setup job with env variables  
	mkdir -p $outdir
	bsub -env "all, INFILE=$j, OUTDIR=$outdir" < mm10_feature_generation_from_sequence.lsf
done


