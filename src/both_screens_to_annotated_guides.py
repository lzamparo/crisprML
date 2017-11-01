import pandas
import os
from gRNA_data_gene_check import load_pickle

# load both data sets
os.chdir(os.path.expanduser("~/projects/crisprML/data/"))

hg_table = pandas.read_csv("raw_features_computed_Cas9_sequences_vs_hg19.csv")
mm_table = pandas.read_csv("raw_features_computed_Cas9_sequences_vs_mm10.csv")

# cut down levenstein dist, summed stuff
hg_table = hg_table[['sequence','Specificity_Score','Occurrences_at_Hamming_0','Occurrences_at_Hamming_1', 'Occurrences_at_Hamming_2','Occurrences_at_Hamming_3','Sum_Hamming_Neighbors']]
mm_table = mm_table[['sequence','Specificity_Score','Occurrences_at_Hamming_0','Occurrences_at_Hamming_1', 'Occurrences_at_Hamming_2','Occurrences_at_Hamming_3','Sum_Hamming_Neighbors']]

# merge tables
merged_scores = pandas.merge(hg_table, mm_table, on="sequence", suffixes=('_hg',"_mm"))

# associate exon targeted with each guide
gRNA_target_coordinates_annotation_dict = load_pickle("gRNA_target_coordinates_annotation_dict.pkl")

# split (on targeted exon), select top 4 (if possible) min(sum_hamming_mm + sum_hamming_hg)??



