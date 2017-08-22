#!/bin/bash

set -e

script="/cbio/cllab/home/aperez/Project/Ventura/ClassTask/code/guideclass/code/gRNA_from_annotations.py"
annotations="/cbio/cllab/home/aperez/Reference_Genomes/mm10/annotation/mm10_annotation_cds_exon_number.bed"
#infile="/cbio/cllab/projects/perez/crispr/all_kmers/mm10/Cas9/mm10_all_kmers.txt"
#infile="/cbio/cllab/projects/perez/crispr/all_kmers/mm10/Cas9/GuideScan_gRNAs.txt"
#infile="/cbio/cllab/projects/perez/crispr/all_kmers/mm10/Cpf1/Cpf1_GuideScan_gRNAs.txt"
infile="/cbio/cllab/projects/perez/crispr/all_kmers/mm10/Cpf1/cpf1_mm10_kmers.txt"
#outdir="/cbio/cllab/home/aperez/Project/Ventura/ClassTask/data/gRNAs_for_screen/mm10/gRNAs_in_mouse_exons/Cas9/GuideScan_database"
#outdir="/cbio/cllab/home/aperez/Project/Ventura/ClassTask/data/gRNAs_for_screen/mm10/gRNAs_in_mouse_exons/Cpf1/GuideScan_database"
outdir="/cbio/cllab/home/aperez/Project/Ventura/ClassTask/data/gRNAs_for_screen/mm10/gRNAs_in_mouse_exons/Cpf1/augmented_annotation"
cpf1=1
output1="gRNA_annotation"
echo "python $script -a $annotations -i $infile -o $outdir --cpf1 $cpf1" #| qsub -l nodes=1:ppn=1,mem=100gb,walltime=72:00:00 -N $output1
