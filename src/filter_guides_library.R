### Characterize a given set of exons wrt how many guides they have per gene, of what quality, etc.
library(data.table)
library(ggplot2)
library(dplyr)


### Get the input for guides vs mm10 and hg19
setwd("/Users/zamparol/projects/crisprML/data")
guide_features_mm10 = data.table(read.csv("guides/raw_features_computed_Cas9_sequences_vs_mm10.csv"))
guide_features_mm10 = guide_features_mm10[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_hg19 = data.table(read.csv("guides/raw_features_computed_Cas9_sequences_hg19.csv"))
guide_features_hg19 = guide_features_hg19[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]

### Get the list of exons that actually overlap guide regions
coding_exons = data.table(read.delim("coding_exons/mm10_first_four_coding_exons_ensembl.bed", sep="\t", header=FALSE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript","gene", "transcript_start", "transcript_end")

### Maybe here is the part where I need to collapse overlapping coding regions?

### Get the map from guides to genes
guides_to_targets = data.table(read.csv("guides/guide_to_target_region.csv", header=FALSE))
colnames(guides_to_targets) = c("guide", "chrom", "start", "end")

### Join the guides based on sequence
feature_tables = merge(x=guide_features_hg19, y=guide_features_mm10, by.x="sequence", by.y="sequence", suffixes=c(".hg", ".mm"))

### Filter guides based on having 0 occurrences in Hg at hamming_0, 1, 2 && 1 occurrrence in MM at hamming_0 but none at hamming_1, 2
filtered_guides = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_1.hg == 0 & Occurrences_at_Hamming_2.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### tidy up environment
rm(guide_features_mm10, guide_features_hg19, feature_tables)

### Join filtered guides targeted gene, transcript, exon, retain info re: human and mouse
setkey(guides_to_targets, chrom, start, end)
setkey(coding_exons, chrom, start, end)
guides_to_exons = foverlaps(guides_to_targets, coding_exons, type="within", nomatch=0L)


### Keep only those guides that are part of the filtered guides set
setkey(filtered_guides, sequence)
setkey(guides_to_exons, guide)
filtered_guides_and_targets = merge(filtered_guides, guides_to_exons, by.x=c("sequence"), by.y=c("guide"))

### What is going on??  How can I have so many guides per gene, exon, exon number??
filtered_guides_and_targets[, num_guides := .N, by = .(gene, transcript, exon)]
ufgs = unique(filtered_guides_and_targets[ , .(num_guides), by = .(gene, transcript, exon)])

### How many (gene, transcript, exon) records do not have guides?
setkey(coding_exons, gene, transcript, exon)
setkey(ufgs, gene, transcript, exon)
merged_w_guides = merge(ufgs, coding_exons, all.y=TRUE)

### How many (gene, transcript, exon) records have guides?
ggplot(ufgs, aes(x=num_guides)) + 
  geom_bar() + 
  xlim(0,50) + 
  ggtitle("Number of high quality guides per gene,transcript,exon") + 
  xlab("Number of guides") + 
  ylab("count") + 
  annotate("text", x = 40, y = 10000, label = paste("(g,t,e) hit:", ufgs[,.N], sep=" ")) + 
  annotate("text", x = 41, y = 9400, label = paste("(g,t,e) not hit:", merged_w_guides[is.na(num_guides), .N], sep=" "))



### Now, for those exons not currently hit by any guides, let's see what are available to us
remaining_exons = merged_w_guides[is.na(num_guides), .(gene, transcript, exon_number, chrom, start, end, total_exons.y)]
setkey(remaining_exons, chrom, start, end)
guides_to_remaining_exons = foverlaps(guides_to_targets, remaining_exons, type="within", nomatch=0L)

### Start winnowing away:
### Need to make sure they are unique in mm10 at Hamming 0, and no hits in hg19 at Hamming 0
setkey(guides_to_remaining_exons, guide)
setkey(feature_tables, sequence)
guides_to_remaining_exons = merge(guides_to_remaining_exons, feature_tables, by.x=c("guide"), by.y=c("sequence"))

filtered_guides_to_remaining_exons = guides_to_remaining_exons[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_0.mm == 1,] 
filtered_guides_to_remaining_exons[, score_summary := Occurrences_at_Hamming_1.hg + Occurrences_at_Hamming_2.hg + Occurrences_at_Hamming_1.mm + Occurrences_at_Hamming_2.mm]

### Group by gene, transcript, exon_number, select the guide with minimal hg_score_summary + mm_score_summary
grouped_filtered_guides_to_remaining_dt = filtered_guides_to_remaining_exons[, .SD[which.min(score_summary)], by=.(gene, transcript, exon_number)]

grouped_filtered_guides_to_remaining <- filtered_guides_to_remaining_exons %>% group_by(gene, transcript, exon_number) %>% top_n(-1, score_summary)
grouped_filtered_guides_to_remaining = as.data.table(ungroup(grouped_filtered_guides_to_remaining))

