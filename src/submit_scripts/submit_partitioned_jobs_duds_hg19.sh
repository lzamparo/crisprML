#! /bin/bash

# submit partitioned specificity score jobs
indir=$1

jobs=$(cat empties.txt)

for j in $jobs
do
	
	# make output directory
	
        outdir="$(dirname $indir)/specificity_output/$j"
	suffix=$(echo $j | cut -f2 -d_)
	infile=$(echo "$indir/Cas9_sequences_fitting_annotation_$suffix.txt")
	# setup job with env variables  
	mkdir -p $outdir
	bsub -env "all, INFILE=$infile, OUTDIR=$outdir" < feature_generation_from_sequence.lsf
done


