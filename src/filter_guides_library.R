### Characterize a given set of exons wrt how many guides they have per gene, of what quality, etc.
library(data.table)
library(ggplot2)

### Get the input for guides vs mm10 and hg19
setwd("/Users/zamparol/projects/crisprML/data")
guide_features_mm10 = data.table(read.csv("raw_features_computed_Cas9_sequences_vs_mm10.csv"))
guide_features_mm10 = guide_features_mm10[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_hg19 = data.table(read.csv("raw_features_computed_Cas9_sequences_vs_hg19.csv"))
guide_features_hg19 = guide_features_hg19[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]

### Get the list of exons that actually overlap guide regions
coding_exons = data.table(read.delim("mm10_annotation_cds_first_four_exons_overlapped_w_coding.bed", sep="\t", header=FALSE))
colnames(coding_exons) = c("chrom","start","end","gene","exon_number", "bah", "strand", "num_overlaps")

### Get the map from guides to genes
guides_to_exons = data.table(read.csv("guide_to_region.csv"))
