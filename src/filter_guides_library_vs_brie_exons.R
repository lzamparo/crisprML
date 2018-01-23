### Characterize a given set of exons wrt how many guides they have per gene, of what quality, etc.
library(data.table)
library(ggplot2)
library(stringr)
library(gridExtra)

### Get the input for guides vs mm10 and hg19
setwd("/Users/zamparol/projects/crisprML/data")
guide_features_mm10 = data.table(fread("guides/main_run/raw_features_computed_Cas9_sequences_vs_mm10.csv"))
guide_features_mm10 = guide_features_mm10[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_mm10 = unique(guide_features_mm10)

guide_features_hg19 = data.table(fread("guides/main_run/raw_features_computed_Cas9_sequences_vs_hg19.csv"))
guide_features_hg19 = guide_features_hg19[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_hg19 = unique(guide_features_hg19)

### Get the list of exons that actually overlap guide regions
coding_exons = data.table(fread("coding_exons/mm10_first_four_filtered_brie_exons.bed", sep="\t", header=TRUE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript","gene", "transcript_start", "transcript_end")
coding_exons = unique(coding_exons)

### Get the map from guides to genes
guides_to_targets = data.table(fread("guides/main_run/guide_to_target_region.csv", header=FALSE))
colnames(guides_to_targets) = c("guide", "chrom", "start", "end")

### Join the guides based on sequence
setkey(guide_features_mm10, sequence)
setkey(guide_features_hg19, sequence)
feature_tables = merge(x=guide_features_hg19, y=guide_features_mm10, by.x="sequence", by.y="sequence", suffixes=c(".hg", ".mm"))

### Filter guides based on having 0 occurrences in Hg at hamming_0, 1, 2 && 1 occurrrence in MM at hamming_0 but none at hamming_1, 2
#filtered_guides = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_1.hg == 0 & Occurrences_at_Hamming_2.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### Filter guides based on having 0 ocurrences in hg19 at Hamming_0 && 1 occurrrence in mm10 at hamming_0 but none at hamming_1, 2
filtered_guides = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

### Eliminate guides which have subsequences that create matches with the restiction enzymes we use
BamHI_filter = "^[G]?ATCC"
BlpI_filter = "^CT[ACGT]{1}AGC"
BstXI_filter = "CCA[ACGT]{6}TG$"
filtered_guides = filtered_guides[(!grepl(BamHI_filter, sequence)) & (!grepl(BlpI_filter, sequence)) & !(grepl(BstXI_filter, sequence)),]

### tidy up environment
rm(guide_features_mm10, guide_features_hg19)

### Join filtered guides targeted gene, transcript, exon, retain info re: human and mouse
setkey(guides_to_targets, chrom, start, end)
setkey(coding_exons, chrom, start, end)
guides_to_exons = foverlaps(guides_to_targets, coding_exons, type="within", nomatch=0L)

### Keep only those guides that are part of the filtered guides set
setkey(filtered_guides, sequence)
setkey(guides_to_exons, guide)
filtered_guides_and_targets = merge(filtered_guides, guides_to_exons, by.x=c("sequence"), by.y=c("guide"))
setnames(filtered_guides_and_targets, c("i.start", "i.end"), c("guide_start", "guide_end"))

### How many guides target each gene?  How many exons per gene?
filtered_guides_and_targets[, num_guides := uniqueN(sequence), by = .(gene)]
genes_and_guides = unique(filtered_guides_and_targets[,.(gene, num_guides)])
covered_genes = genes_and_guides[num_guides >= 4, gene]
filtered_guides_and_targets[, exon_count := uniqueN(exon), by = gene]

### What about just the guide which targets each exon, but with minimal Hamming_3.mm + Hamming_3.hg?
filtered_guides_and_targets[, hamming_3_sum := Occurrences_at_Hamming_3.mm + Occurrences_at_Hamming_3.hg]

### Break down filtered_guides_and_targets by gene with 4,3,2,1 exons
four_exons_fgt = filtered_guides_and_targets[gene %in% covered_genes & exon_count >= 4,]
three_exons_fgt = filtered_guides_and_targets[gene %in% covered_genes & exon_count == 3,]
two_exons_fgt = filtered_guides_and_targets[gene %in% covered_genes & exon_count == 2,]
one_exons_fgt = filtered_guides_and_targets[gene %in% covered_genes & exon_count == 1,]


### Four: choose one guide / exon, each with minimum hamming_3_sum
four_exon_fgt_guides = four_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
four_exon_fgt_guides = four_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, hamming_3_sum, guide_start, guide_end)]

### Three: choose one guide / exon, and then the minimum hamming_3_sum of remaining guides
three_exon_fgt_guides = three_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
chosen_guides = three_exon_fgt_guides[,sequence]
additionals = three_exons_fgt[!(sequence %in% chosen_guides), .SD[head(order(hamming_3_sum),1)], by = .(transcript)]
three_exon_fgt_guides = data.table(rbind(three_exon_fgt_guides,additionals))
three_exon_fgt_guides = three_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, hamming_3_sum, guide_start, guide_end)]

### Two: choose one guide / exon, and then the smallest two hamming_3_sum of remaining guides
two_exon_fgt_guides = two_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
chosen_guides = two_exon_fgt_guides[,sequence]
additionals = two_exons_fgt[!(sequence %in% chosen_guides), .SD[head(order(hamming_3_sum),2)], by = .(transcript)]
two_exon_fgt_guides = data.table(rbind(two_exon_fgt_guides,additionals))
two_exon_fgt_guides = two_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, hamming_3_sum, guide_start, guide_end)]


### One: choose the four guides with minimal hamming_3_sum
one_exon_fgt_guides = one_exons_fgt[, .SD[head(order(hamming_3_sum),4)], by = .(exon)]
one_exon_fgt_guides = one_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, hamming_3_sum, guide_start, guide_end)]


### join, write out guides
all_fgt_guides = data.table(rbind(four_exon_fgt_guides,three_exon_fgt_guides,two_exon_fgt_guides,one_exon_fgt_guides))
write.csv(all_fgt_guides,file = "guides/main_run/all_fgt_guides_nobaddies.csv", row.names=FALSE)

### write out exons needing more guides
all_fgt_guides[, guides_per_tx := .N, by=transcript]
txs_per_gene = unique(all_fgt_guides[,.(transcript, guides_per_tx)])
all_gt = unique(coding_exons[, gene, by=transcript])

guides_per_tx = merge(all_gt, txs_per_gene, by=c("transcript"), all.x=TRUE)
guides_per_tx[is.na(guides_per_tx), guides_per_tx := 0]
write.csv(guides_per_tx[guides_per_tx < 4,], "coding_exons/txs_needing_more_guides.csv", row.names=FALSE)

### Plots which present the number of guides / gene gained by relaxing constaints on hg19
### N.B: do not run in sequence. The strict selection results and lax seletion results are derived from separate runs.
setkey(genes_and_guides, gene)
all_gt = unique(coding_exons[, gene, by=transcript])
setkey(all_gt, gene)
all_gt_merged = merge(all_gt, genes_and_guides, by=c("gene"), all.x=TRUE)
all_gt_merged[is.na(num_guides), num_guides := 0]

p1 = ggplot(all_gt_merged, aes(x=num_guides)) + geom_bar(position="dodge") + xlab("# guides per gene") + ylab("Number of genes") + coord_cartesian(xlim = c(0,20)) +  ggtitle("Count of number of guides per gene (zoom: [0,20])")

# selected guides, strict criteria
strict_selection_results = all_fgt_guides[, .N, by=transcript]
strict_selection_results[, selection := "strict"]

# selected guides, lax criteria
lax_selection_results = all_fgt_guides[, .N, by=transcript]
lax_selection_results[, selection := "relaxed"]

guide_improvement = data.table(rbind(strict_selection_results, lax_selection_results))

p2 = ggplot(guide_improvement[N == 4,], aes(x=N, fill=selection)) + geom_bar(position = "dodge") + xlab("Guide thresholds") +
  ylab("Number of genes") + 
  coord_cartesian(ylim = c(17000,19500)) + 
  ggtitle("Number of genes with four guides") + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())



