#!/bin/bash
#BSUB -J mm10_guide_test
#BSUB -n 1
#BSUB -R rusage[mem=100]
#BSUB -W 4:00
#BSUB -o %J.stdout
#BSUB -eo %J.stderr

# set up env, cd to source
annotations="/data/leslie/zamparol/crisprML/genomes/mm10/mm10_annotation_cds_exon_number.bed"
infile="/data/leslie/zamparol/crisprML/kmers/mm10/Cas9/mm10_all_kmers.txt"
outdir="/data/leslie/zamparol/crisprML/results/gRNAs/mm10/"


cd ~/projects/crisprML/src
python gRNA_from_annotations.py -a $annotations -i $infile -o $outdir 
