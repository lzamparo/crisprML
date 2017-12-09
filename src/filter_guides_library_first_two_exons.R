### Characterize a given set of exons wrt how many guides they have per gene, of what quality, etc.
library(data.table)
library(ggplot2)
library(dplyr)


### Get the input for guides vs mm10 and hg19
setwd("/Users/zamparol/projects/crisprML/data")
guide_features_mm10 = data.table(fread("guides/raw_features_computed_Cas9_sequences_vs_mm10.csv"))
guide_features_mm10 = guide_features_mm10[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_hg19 = data.table(fread("guides/raw_features_computed_Cas9_sequences_vs_hg19.csv"))
guide_features_hg19 = guide_features_hg19[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]

### Get the list of exons that actually overlap guide regions
coding_exons = data.table(fread("coding_exons/mm10_first_four_coding_exons_ensembl.bed", sep="\t", header=TRUE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript","gene", "transcript_start", "transcript_end")

### Select the first two coding exons from each trancripts
first_two_coding_exons = coding_exons[, .SD[order(start)][1:2,], by=.(transcript)]

### the above works, but misses out on those transcripts which are mono-exonic
mono_exonic_txs = coding_exons[transcript %in% first_two_coding_exons[is.na(chrom) & is.na(start) & is.na(end),transcript], ]
first_two_coding_exons = first_two_coding_exons[!is.na(chrom),]
first_two_coding_exons = data.table(rbind(first_two_coding_exons,mono_exonic_txs))

### annotate how many exons are shared by > 1 tx
first_two_coding_exons[,exon_count := .N,by = exon]
# first_two_coding_exons[,.N] - first_two_coding_exons[exon_count > 1, .N]
#[1] 131518

### Get the map from guides to genes
guides_to_targets = data.table(fread("guides/guide_to_target_region.csv", header=FALSE))
colnames(guides_to_targets) = c("guide", "chrom", "start", "end")

### Join the guides based on sequence
feature_tables = merge(x=guide_features_hg19, y=guide_features_mm10, by.x="sequence", by.y="sequence", suffixes=c(".hg", ".mm"))

### Filter guides based on having 0 occurrences in Hg at hamming_0, 1, 2 && 1 occurrrence in MM at hamming_0 but none at hamming_1, 2
filtered_guides = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_1.hg == 0 & Occurrences_at_Hamming_2.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### tidy up environment
rm(guide_features_mm10, guide_features_hg19, coding_exons)

### Join filtered guides targeted gene, transcript, exon, retain info re: human and mouse
setkey(guides_to_targets, chrom, start, end)
setkey(first_two_coding_exons, chrom, start, end)
guides_to_exons = foverlaps(guides_to_targets, first_two_coding_exons, type="within", nomatch=0L)

### Keep only those guides that are part of the filtered guides set
setkey(filtered_guides, sequence)
setkey(guides_to_exons, guide)
filtered_guides_and_targets = merge(filtered_guides, guides_to_exons, by.x=c("sequence"), by.y=c("guide"))
setnames(filtered_guides_and_targets, c("i.start", "i.end"), c("guide_start", "guide_end"))
#saveRDS(filtered_guides_and_targets, file="../results/library/more_specific_guides_and_targets.rds")

### How many guides target each exon?
filtered_guides_and_targets[, num_guides := .N, by = .(exon)]

### What about just the guide which targets each exon, but with minimal Hamming_3 in mm?
zug = filtered_guides_and_targets[, .SD[which.min(Occurrences_at_Hamming_3.mm)], by = .(exon)]

