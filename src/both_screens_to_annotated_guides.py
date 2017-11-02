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

# select those guides with no human hits for 0,1,2 mismatches && no mouse hits for 1,2 mismatches
guide_superset = merged_scores[(merged_scores['Occurrences_at_Hamming_0_hg'] == 0.0) & (merged_scores['Occurrences_at_Hamming_1_hg'] == 0.0) & (merged_scores['Occurrences_at_Hamming_2_hg'] == 0.0) & (merged_scores['Occurrences_at_Hamming_0_mm'] == 1.0) & (merged_scores['Occurrences_at_Hamming_1_mm'] == 0.0) & (merged_scores['Occurrences_at_Hamming_2_mm'] == 0.0)]

# associate exon targeted with each guide
gRNA_target_coordinates_annotation_dict = load_pickle("gRNA_target_coordinates_annotation_dict.pkl")


#["chr2:35084609-35084609_['AI182371']"]
def to_csv_rec(elem):
    chrom = elem.split(":")[0]
    rest = elem.split(":")[1]
    start = rest.split("-")[0]
    raggy = rest.split("-")[1]
    end = raggy.split("_")[0]
    return ','.join([chrom,start,end])


    


