#!/bin/bash

set -e

# dataset sgRNA fields

script="/cbio/cllab/home/aperez/Project/Ventura/ClassTask/code/specificity_score_distance_neighbors.py"
#infile="/cbio/cllab/projects/perez/classtask/data/gRNAs_for_screen/NGG_filtered_gRNAs/SacCerv_all_gRNAs_NGG.txt"
outdir="/cbio/cllab/home/aperez/Project/Ventura/ClassTask/data/gRNAs_for_screen/mm10/gRNAs_in_mouse_exons/Cpf1/augmented_annotation/split/human_trie_processed"
kmers="/cbio/cllab/projects/perez/crispr/kmers_counted/Cpf1/hg38/hg38_all_cpf1_kmers_counted.txt"
trie="/cbio/cllab/home/aperez/cpf1_hg38/cpf1_hg38_trie.dat"
#kmers="/cbio/cllab/projects/perez/crispr/kmers_counted/Cpf1/mm10/mm10_all_cpf1_kmers_counted.txt"
#trie="/cbio/cllab/home/aperez/cpf1_mm10/cpf1_mm10_trie.dat"
#kmers="/cbio/cllab/projects/perez/crispr/kmers_counted/mm10/mm10_all_kmers_counted.txt" #mouse
#trie="/cbio/cllab/projects/pritykin/CRISPR/mm10/mm10_all/mm10_all_trie.dat" #mouse
#kmers="/cbio/cllab/projects/perez/crispr/kmers_counted/hg38/hg38_all_kmers_counted.txt" #human
#trie="/cbio/cllab/projects/pritykin/CRISPR/hg38/hg38_all/hg38_all_trie.dat" # human
mismatch="/cbio/cllab/home/aperez/Scripts/Project_Scripts/Ventura_CRISPR/active/crispr-project/guidescan-crispr/guidescan/CFD_scoring/mismatch_score.pkl"
pam="/cbio/cllab/home/aperez/Scripts/Project_Scripts/Ventura_CRISPR/active/crispr-project/guidescan-crispr/guidescan/CFD_scoring/pam_scores.pkl"
header="True"
sequence_field=0 #default 0
cpf1=0
basedir="/cbio/cllab/home/aperez/Project/Ventura/ClassTask/data/gRNAs_for_screen/mm10/gRNAs_in_mouse_exons/Cpf1/augmented_annotation/split"
export PYTHONPATH=/cbio/grlab/share/software/lib/python2.7/site-packages:$PYTHONPATH

output1="featurize_from_sequence"

for FILE in $(find $basedir -name cpf1_augmented_new_IT*)
	do
		echo "python $script -i $FILE -o $outdir -k $kmers -t $trie -m $mismatch -p $pam --sequence_field $sequence_field --cpf1 $cpf1" #| qsub -l nodes=1:ppn=1,mem=200gb,walltime=75:00:00 -N $output1
	done
