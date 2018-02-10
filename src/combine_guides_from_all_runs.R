library(data.table)
library(ggplot2)

# load the main run guides and target maps
setwd("/Users/zamparol/projects/crisprML/data")
### Get the list of exons that actually overlap guide regions
coding_exons = data.table(fread("coding_exons/mm10_first_four_filtered_brie_exons.bed", sep="\t", header=TRUE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript","gene", "transcript_start", "transcript_end")
coding_exons = unique(coding_exons)

##### Main run
guide_features_mm10 = data.table(fread("guides/main_run/raw_features_computed_Cas9_sequences_vs_mm10.csv"))
guide_features_mm10 = guide_features_mm10[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_mm10 = unique(guide_features_mm10)

guide_features_hg19 = data.table(fread("guides/main_run/raw_features_computed_Cas9_sequences_vs_hg19.csv"))
guide_features_hg19 = guide_features_hg19[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_hg19 = unique(guide_features_hg19)

### Get the map from guides to genes
guides_to_targets = data.table(fread("guides/main_run/guide_to_target_region.csv", header=FALSE))
colnames(guides_to_targets) = c("guide", "chrom", "start", "end")

### Join the guides based on sequence
setkey(guide_features_mm10, sequence)
setkey(guide_features_hg19, sequence)
feature_tables = merge(x=guide_features_hg19, y=guide_features_mm10, by.x="sequence", by.y="sequence", suffixes=c(".hg", ".mm"))

### Filter guides based on having 0 ocurrences in hg19 at Hamming_0 && 1 occurrrence in mm10 at hamming_0 but none at hamming_1, 2
filtered_guides = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0,]

### Eliminate guides which have subsequences that create matches with the restiction enzymes we use
BamHI_filter_emergent = "^[G]?ATCC"
BlpI_filter_emergent = "^CT[ACGT]{1}AGC"
BstXI_filter_emergent = "CCA[ACGT]{6}TG$"

BstXI_filter_exact = "CCA[ACGT]{6}TGG"
BamHI_filter_exact = "GGATCC"
BlpI_filter_exact = "GCT[ACGT]{1}AGC"

filtered_guides = filtered_guides[(str_count(sequence, BamHI_filter_exact) == 0) & (str_count(sequence, BlpI_filter_exact) == 0) & (str_count(sequence, BstXI_filter_exact) == 0) & (str_count(sequence, BamHI_filter_emergent) == 0) & (str_count(sequence, BlpI_filter_emergent) == 0) & (str_count(sequence, BstXI_filter_emergent) == 0),]

### tidy up environment
rm(guide_features_mm10, guide_features_hg19)

### Get the map from guides to genes
guides_to_targets = data.table(fread("guides/main_run/guide_to_target_region.csv", header=FALSE))
colnames(guides_to_targets) = c("guide", "chrom", "start", "end")

setkey(guides_to_targets, chrom, start, end)
setkey(coding_exons, chrom, start, end)
guides_to_exons = foverlaps(guides_to_targets, coding_exons, type="within", nomatch=0L)

### Keep only those guides that are part of the filtered guides set
setkey(filtered_guides, sequence)
setkey(guides_to_exons, guide)
filtered_guides_and_targets_main = merge(filtered_guides, guides_to_exons, by.x=c("sequence"), by.y=c("guide"))
setnames(filtered_guides_and_targets_main, c("i.start", "i.end"), c("guide_start", "guide_end"))

##### Low coverage run

guide_features_mm10 = data.table(fread("guides/later_needy_exon_run/raw_features_computed_Cas9_sequences_vs_mm10.csv"))
guide_features_mm10 = guide_features_mm10[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_mm10 = unique(guide_features_mm10)

guide_features_hg19 = data.table(fread("guides/later_needy_exon_run/raw_features_computed_Cas9_sequences_vs_hg19.csv"))
guide_features_hg19 = guide_features_hg19[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_hg19 = unique(guide_features_hg19)

guides_to_targets = data.table(fread("guides/later_needy_exon_run/guide_to_target_region.csv", header=FALSE))
colnames(guides_to_targets) = c("guide", "chrom", "start", "end","strand")

### Join the guides based on sequence
setkey(guide_features_mm10, sequence)
setkey(guide_features_hg19, sequence)
feature_tables = merge(x=guide_features_hg19, y=guide_features_mm10, by.x="sequence", by.y="sequence", suffixes=c(".hg", ".mm"))

filtered_guides_1 = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_1.hg == 0 & Occurrences_at_Hamming_2.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### Filter guides based on having 0 ocurrences in hg19 at Hamming_0 && 1 occurrrence in mm10 at hamming_0 but none at hamming_1, 2
filtered_guides_2 = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### Filter guides based only on restrictions in mm10
filtered_guides_3 = feature_tables[Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### Filter guides based on relaxed restrictions in mm10
filtered_guides_4 = feature_tables[Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0,]

### Eliminate guides which have subsequences that create matches with the restiction enzymes we use
filtered_guides = unique(data.table(rbind(filtered_guides_1,filtered_guides_2,filtered_guides_3,filtered_guides_4)))
filtered_guides = filtered_guides[(str_count(sequence, BamHI_filter_exact) == 0) & (str_count(sequence, BlpI_filter_exact) == 0) & (str_count(sequence, BstXI_filter_exact) == 0) & (str_count(sequence, BamHI_filter_emergent) == 0) & (str_count(sequence, BlpI_filter_emergent) == 0) & (str_count(sequence, BstXI_filter_emergent) == 0),]

### tidy up environment
rm(guide_features_mm10, guide_features_hg19)

setkey(guides_to_targets, chrom, start, end)
setkey(coding_exons, chrom, start, end)
guides_to_exons = foverlaps(guides_to_targets, coding_exons, type="within", nomatch=0L)

### Start by assigning only those guides that are part of the filtered guides set
setkey(filtered_guides, sequence)
setkey(guides_to_exons, guide)
filtered_guides_and_targets_needy = merge(filtered_guides, guides_to_exons, by.x=c("sequence"), by.y=c("guide"))
setnames(filtered_guides_and_targets_needy, c("i.start", "i.end"), c("guide_start", "guide_end"))

##### missing exon run

guide_features_mm10 = data.table(fread("guides/missing_exons_run/raw_features_computed_Cas9_sequences_vs_mm10.csv"))
guide_features_mm10 = guide_features_mm10[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_mm10 = unique(guide_features_mm10)

guide_features_hg19 = data.table(fread("guides/missing_exons_run/raw_features_computed_Cas9_sequences_vs_hg19.csv"))
guide_features_hg19 = guide_features_hg19[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_hg19 = unique(guide_features_hg19)

guides_to_targets = data.table(fread("guides/missing_exons_run/guide_to_target_region.csv", header=FALSE))
colnames(guides_to_targets) = c("guide", "chrom", "start", "end","strand")

### Join the guides based on sequence
setkey(guide_features_mm10, sequence)
setkey(guide_features_hg19, sequence)
feature_tables = merge(x=guide_features_hg19, y=guide_features_mm10, by.x="sequence", by.y="sequence", suffixes=c(".hg", ".mm"))

### Filter guides based on having 0 occurrences in Hg at hamming_0, 1, 2 && 1 occurrrence in MM at hamming_0 but none at hamming_1, 2
filtered_guides_1 = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_1.hg == 0 & Occurrences_at_Hamming_2.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### Filter guides based on having 0 ocurrences in hg19 at Hamming_0 && 1 occurrrence in mm10 at hamming_0 but none at hamming_1, 2
filtered_guides_2 = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### Filter guides based only on restrictions in mm10
filtered_guides_3 = feature_tables[Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### Filter guides based on relaxed restrictions in mm10
filtered_guides_4 = feature_tables[Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0,]

### Eliminate guides which have subsequences that create matches with the restiction enzymes we use
filtered_guides = unique(data.table(rbind(filtered_guides_1,filtered_guides_2,filtered_guides_3,filtered_guides_4)))
filtered_guides = filtered_guides[(str_count(sequence, BamHI_filter_exact) == 0) & (str_count(sequence, BlpI_filter_exact) == 0) & (str_count(sequence, BstXI_filter_exact) == 0) & (str_count(sequence, BamHI_filter_emergent) == 0) & (str_count(sequence, BlpI_filter_emergent) == 0) & (str_count(sequence, BstXI_filter_emergent) == 0),]

### tidy up environment
rm(guide_features_mm10, guide_features_hg19)

setkey(guides_to_targets, chrom, start, end)
setkey(coding_exons, chrom, start, end)
guides_to_exons = foverlaps(guides_to_targets, coding_exons, type="within", nomatch=0L)

### Start by assigning only those guides that are part of the filtered guides set
setkey(filtered_guides, sequence)
setkey(guides_to_exons, guide)
filtered_guides_and_targets_missing = merge(filtered_guides, guides_to_exons, by.x=c("sequence"), by.y=c("guide"))
setnames(filtered_guides_and_targets_missing, c("i.start", "i.end"), c("guide_start", "guide_end"))

##### merge all filtered guides and targets, match up any guides for low coverage genes
filtered_guides_and_targets_needy = filtered_guides_and_targets_needy[,.(sequence, Specificity_Score.hg, Occurrences_at_Hamming_0.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_3.hg, Specificity_Score.mm, Occurrences_at_Hamming_0.mm, Occurrences_at_Hamming_1.mm, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_3.mm, chrom, start, end, strand, exon, transcript, gene, transcript_start, transcript_end, guide_start, guide_end)]
filtered_guides_and_targets_missing = filtered_guides_and_targets_missing[,.(sequence, Specificity_Score.hg, Occurrences_at_Hamming_0.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_3.hg, Specificity_Score.mm, Occurrences_at_Hamming_0.mm, Occurrences_at_Hamming_1.mm, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_3.mm, chrom, start, end, strand, exon, transcript, gene, transcript_start, transcript_end, guide_start, guide_end)]
all_filtered_guides_and_targets = unique(data.table(rbind(filtered_guides_and_targets_main, filtered_guides_and_targets_needy, filtered_guides_and_targets_missing)))

### load the low coverage genes
low_genes = data.table(fread("fewer_than_four_genes.csv"))
all_against_low_genes = all_filtered_guides_and_targets[exon %in% low_genes[,exon] & transcript %in% low_genes[,transcript],]

##### grab guides against genes needing more representaiton among 

all_against_low_genes[, exon_count := uniqueN(exon), by = gene]
all_against_low_genes[, hamming_3_sum := Occurrences_at_Hamming_3.mm + Occurrences_at_Hamming_3.hg]

### Break down all_against_low_genes by gene with 4,3,2,1 exons
four_exons_fgt = all_against_low_genes[exon_count >= 4,]
three_exons_fgt = all_against_low_genes[exon_count == 3,]
two_exons_fgt = all_against_low_genes[exon_count == 2,]
one_exons_fgt = all_against_low_genes[exon_count == 1,]

### Four: choose one guide / exon, each with minimum hamming_3_sum, first four exons by position
four_exon_fgt_guides = four_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
four_exon_fgt_guides = four_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]
four_exon_fgt_guides = four_exon_fgt_guides[,.SD[head(order(start),4)], by = .(transcript)]

### Three: choose one guide / exon, and then the minimum hamming_3_sum of remaining guides
three_exon_fgt_guides = three_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
chosen_guides = three_exon_fgt_guides[,sequence]
additionals = three_exons_fgt[!(sequence %in% chosen_guides), .SD[head(order(hamming_3_sum),1)], by = .(transcript)]
three_exon_fgt_guides = data.table(rbind(three_exon_fgt_guides,additionals))
three_exon_fgt_guides = three_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]

### Two: choose one guide / exon, and then the smallest two hamming_3_sum of remaining guides
two_exon_fgt_guides = two_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
chosen_guides = two_exon_fgt_guides[,sequence]
additionals = two_exons_fgt[!(sequence %in% chosen_guides), .SD[head(order(hamming_3_sum),2)], by = .(transcript)]
two_exon_fgt_guides = data.table(rbind(two_exon_fgt_guides,additionals))
two_exon_fgt_guides = two_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]


### One: choose the four guides with minimal hamming_3_sum
one_exon_fgt_guides = one_exons_fgt[, .SD[head(order(hamming_3_sum),4)], by = .(exon)]
one_exon_fgt_guides = one_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]

### combine the data from all runs 
all_fgt_guides = data.table(rbind(four_exon_fgt_guides,three_exon_fgt_guides,two_exon_fgt_guides,one_exon_fgt_guides))

### Make really sure that no guide has any of the restriction enzyme matches (or partial matches).  
### N.B: turns out it was a subtle difference between str_count and grepl
all_fgt_guides[str_count(sequence, BstXI_filter_exact) > 0,.N]
all_fgt_guides[str_count(sequence, BamHI_filter_exact) > 0,.N]
all_fgt_guides[str_count(sequence, BlpI_filter_exact) > 0,.N]
all_fgt_guides[str_count(sequence, BamHI_filter_emergent) > 0,.N]
all_fgt_guides[str_count(sequence, BstXI_filter_emergent) > 0,.N]
all_fgt_guides[str_count(sequence, BlpI_filter_emergent) > 0,.N]


setwd("~/projects/crisprML/results/library/")
write.csv(all_fgt_guides, "rescue_guides_for_low_genes.csv", row.names=FALSE)

### how many do we save?
all_fgt_guides[,guides_per_tx := uniqueN(sequence), by=transcript]
p1 = ggplot(all_fgt_guides[,guides_per_tx, by=transcript], aes(x=guides_per_tx)) + geom_bar(position="dodge") + xlab("# guides per gene") + ylab("Number of genes") + coord_cartesian(xlim = c(0,5)) +  ggtitle("Count of number of guides per gene")


### Re-do guide selection for main run guides
filtered_guides_and_targets_main[, exon_count := uniqueN(exon), by = gene]
filtered_guides_and_targets_main[, hamming_3_sum := Occurrences_at_Hamming_3.mm + Occurrences_at_Hamming_3.hg]

### Break down filtered_guides_and_targets_main by gene with 4,3,2,1 exons
four_exons_fgt = filtered_guides_and_targets_main[exon_count >= 4,]
three_exons_fgt = filtered_guides_and_targets_main[exon_count == 3,]
two_exons_fgt = filtered_guides_and_targets_main[exon_count == 2,]
one_exons_fgt = filtered_guides_and_targets_main[exon_count == 1,]

### Four: choose one guide / exon, each with minimum hamming_3_sum, first four exons by position
four_exon_fgt_guides = four_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
four_exon_fgt_guides = four_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]
four_exon_fgt_guides = four_exon_fgt_guides[,.SD[head(order(start),4)], by = .(transcript)]

### Three: choose one guide / exon, and then the minimum hamming_3_sum of remaining guides
three_exon_fgt_guides = three_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
chosen_guides = three_exon_fgt_guides[,sequence]
additionals = three_exons_fgt[!(sequence %in% chosen_guides), .SD[head(order(hamming_3_sum),1)], by = .(transcript)]
three_exon_fgt_guides = data.table(rbind(three_exon_fgt_guides,additionals))
three_exon_fgt_guides = three_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]

### Two: choose one guide / exon, and then the smallest two hamming_3_sum of remaining guides
two_exon_fgt_guides = two_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
chosen_guides = two_exon_fgt_guides[,sequence]
additionals = two_exons_fgt[!(sequence %in% chosen_guides), .SD[head(order(hamming_3_sum),2)], by = .(transcript)]
two_exon_fgt_guides = data.table(rbind(two_exon_fgt_guides,additionals))
two_exon_fgt_guides = two_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]


### One: choose the four guides with minimal hamming_3_sum
one_exon_fgt_guides = one_exons_fgt[, .SD[head(order(hamming_3_sum),4)], by = .(exon)]
one_exon_fgt_guides = one_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]

### combine the data from all runs 
all_fgt_guides = data.table(rbind(four_exon_fgt_guides,three_exon_fgt_guides,two_exon_fgt_guides,one_exon_fgt_guides))

### Test one more time to see if somehow I do not filter out bogus guides
all_fgt_guides[str_count(sequence, BstXI_filter_exact) > 0,.N]
all_fgt_guides[str_count(sequence, BamHI_filter_exact) > 0,.N]
all_fgt_guides[str_count(sequence, BlpI_filter_exact) > 0,.N]
all_fgt_guides[str_count(sequence, BamHI_filter_emergent) > 0,.N]
all_fgt_guides[str_count(sequence, BstXI_filter_emergent) > 0,.N]
all_fgt_guides[str_count(sequence, BlpI_filter_emergent) > 0,.N]

setwd("~/projects/crisprML/results/library/")
write.csv(all_fgt_guides, "main_run_four_guides_per_gene_final.csv", row.names=FALSE)
