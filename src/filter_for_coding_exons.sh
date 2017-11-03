#!/bin/bash

# filter those exons targeted in the screen for only coding exons

# prefix for where data lies
prefix=/data/leslie/zamparol/crisprML/genomes/mm10
echo "data lies below: $prefix"

# Intersect the exons with list of coding sequence spans from UCSC 
#echo "performing intersection..."
#bedtools intersect -f 1.0 -c -a /home/zamparol/projects/crisprML/data/annotations/mm10_exons_by_symbol_fol.bed -b $prefix/mm10_cds_reordered_sorted.tsv > $prefix/mm10_annotation_cds_first_four_exons_overlapped_w_coding.bed

overlapcounts=$prefix/mm10_annotation_cds_first_four_exons_coding_overlap_50.bed

# Keep only those exons which overlap with coding sequence
cat $overlapcounts | egrep "[1-9]+$" > $prefix/mm10_annotation_cds_first_four_exons_overlapped_w_coding.bed
echo "Kept $(wc -l $prefix/mm10_annotation_cds_first_four_exons_overlapped_w_coding.bed) exons that lie in coding regions"
 
