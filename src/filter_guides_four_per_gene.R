### Characterize a given set of exons wrt how many guides they have per gene, of what quality, etc.

library(data.table)
library(ggplot2)
library(dplyr)


### Get the input for guides vs mm10 and hg19
setwd("/Users/zamparol/projects/crisprML/data")
guide_features_mm10 = data.table(fread("guides/raw_features_computed_Cas9_sequences_vs_mm10.csv"))
guide_features_mm10 = guide_features_mm10[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_mm10 = unique(guide_features_mm10)

guide_features_hg19 = data.table(fread("guides/raw_features_computed_Cas9_sequences_vs_hg19.csv"))
guide_features_hg19 = guide_features_hg19[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_hg19 = unique(guide_features_hg19)

### Get the list of exons that actually overlap guide regions
coding_exons = data.table(fread("coding_exons/mm10_first_four_filtered_brie_exons.bed", sep="\t", header=TRUE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript","gene", "transcript_start", "transcript_end")
coding_exons = unique(coding_exons)

### Get the map from guides to genes
guides_to_targets = data.table(fread("guides/guide_to_target_region.csv", header=FALSE))
colnames(guides_to_targets) = c("guide", "chrom", "start", "end")

### Join the guides based on sequence
setkey(guide_features_mm10, sequence)
setkey(guide_features_hg19, sequence)
feature_tables = merge(x=guide_features_hg19, y=guide_features_mm10, by.x="sequence", by.y="sequence", suffixes=c(".hg", ".mm"))

### Filter guides based on having 0 occurrences in Hg at hamming_0, 1, 2 && 1 occurrrence in MM at hamming_0 but none at hamming_1, 2
filtered_guides = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_1.hg == 0 & Occurrences_at_Hamming_2.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### tidy up environment
rm(guide_features_mm10, guide_features_hg19)

### Need our set of guides to satisfy two criteria:
# 1) We need to cover each gene 
# 2) We need four guides per gene

### Join filtered guides targeted gene, transcript, exon, retain info re: human and mouse
setkey(guides_to_targets, chrom, start, end)
setkey(coding_exons, chrom, start, end)
guides_to_exons = foverlaps(guides_to_targets, coding_exons, type="within", nomatch=0L)

### Start by assigning only those guides that are part of the filtered guides set
setkey(filtered_guides, sequence)
setkey(guides_to_exons, guide)
filtered_guides_and_targets = merge(filtered_guides, guides_to_exons, by.x=c("sequence"), by.y=c("guide"))
setnames(filtered_guides_and_targets, c("i.start", "i.end"), c("guide_start", "guide_end"))

### How many genes have at least four guides?
#filtered_guides_and_targets[, num_guides := .N, by = .(exon)]

### Find the genes which need more guides


### Can we assign more quality guides to those exons such that we satisfy (2)?


### If not, let's find all exons for the remaining genes, and run a separate experiment


### Finally, let's relax the constraints on guides for the human genome; how many genes does this allow us to cover by satisfying conditions (1), (2) ?


### Write out all the guides we have chosen



