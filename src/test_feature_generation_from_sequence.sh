#!/bin/bash


# dataset sgRNA fields
input="/data/leslie/zamparol/crisprML/results/gRNAs/mm10/Cas9_sequences_fitting_annotation.txt"
outdir="/data/leslie/zamparol/crisprML/results/gRNAs/mm10/"
kmers="/data/leslie/zamparol/crisprML/kmers/hg38/Cas9/hg38_kmers.db"
trie="/data/leslie/zamparol/crisprML/tries/mm10_all_trie.dat"
mismatch="/data/leslie/zamparol/crisprML/CFD_scoring/mismatch_score.pkl"
pam="/data/leslie/zamparol/crisprML/CFD_scoring/pam_scores.pkl"


cd ~/projects/crisprML/src

set +o nounset
source activate py2
python specificity_score_distance_neighbors.py -i $input -o $outdir -k $kmers -t $trie -m $mismatch -p $pam

