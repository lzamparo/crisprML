### Characterize a given set of exons wrt how many guides they have per gene, of what quality, etc.
library(data.table)
library(ggplot2)


### Get the input for guides vs mm10 and hg19
setwd("/Users/zamparol/projects/crisprML/data")
guide_features_mm10 = data.table(fread("guides/later_needy_exon_run/raw_features_computed_Cas9_sequences_vs_mm10.csv"))
guide_features_mm10 = guide_features_mm10[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_mm10 = unique(guide_features_mm10)

guide_features_hg19 = data.table(fread("guides/later_needy_exon_run/raw_features_computed_Cas9_sequences_vs_hg19.csv"))
guide_features_hg19 = guide_features_hg19[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_hg19 = unique(guide_features_hg19)

### Get the list of exons that actually overlap guide regions
coding_exons = data.table(fread("coding_exons/mm10_needy_exons.bed", sep="\t", header=TRUE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript", "transcript_start", "transcript_end", "gene")
coding_exons = unique(coding_exons)

### Get the map from guides to genes
guides_to_targets = data.table(fread("guides/later_needy_exon_run/guide_to_target_region.csv", header=FALSE))
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
BamHI_filter = "^[G]?ATCC"
BlpI_filter = "^CT[ACGT]{1}AGC"
BstXI_filter = "CCA[ACGT]{6}TG$"

filtered_guides = unique(data.table(rbind(filtered_guides_1,filtered_guides_2,filtered_guides_3,filtered_guides_4)))
filtered_guides = filtered_guides[(!grepl(BamHI_filter, sequence)) & (!grepl(BlpI_filter, sequence)) & !(grepl(BstXI_filter, sequence)),]

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


### How many guides target each gene?  How many exons per gene?
### Also help score each guide by the sum of their 3-mismatch matches in human and mouse
filtered_guides_and_targets[, exon_count := uniqueN(exon), by = gene]
filtered_guides_and_targets[, hamming_3_sum := Occurrences_at_Hamming_3.mm + Occurrences_at_Hamming_3.hg]


### Break down filtered_guides_and_targets by gene with 4,3,2,1 exons
four_exons_fgt = filtered_guides_and_targets[exon_count >= 4,]
three_exons_fgt = filtered_guides_and_targets[exon_count == 3,]
two_exons_fgt = filtered_guides_and_targets[exon_count == 2,]
one_exons_fgt = filtered_guides_and_targets[exon_count == 1,]

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


### Validate number of guides / gene
four_exon_fgt_guides[,guides_per_gene := uniqueN(sequence),by=gene]
four_exon_fgt_guides[guides_per_gene < 4,]

three_exon_fgt_guides[, guides_per_gene := uniqueN(sequence),by=gene]
three_exon_fgt_guides[guides_per_gene < 4, uniqueN(gene)]

two_exon_fgt_guides[,guides_per_gene := uniqueN(sequence),by=gene]
two_exon_fgt_guides[guides_per_gene < 4, uniqueN(gene)]

one_exon_fgt_guides[,guides_per_gene := uniqueN(sequence),by=gene]
one_exon_fgt_guides[guides_per_gene < 4, uniqueN(gene)]

### Figure out how to show discrepancy:



### Combine all guides together, write out guides
all_fgt_guides = data.table(rbind(four_exon_fgt_guides,three_exon_fgt_guides,two_exon_fgt_guides,one_exon_fgt_guides))
#write.csv(all_fgt_guides,file = "guides/needy_exon_run/all_fgt_guides_nobaddies.csv", row.names=FALSE)

### Plots which present the number of guides / gene gained by relaxing constaints on hg19
### N.B: do not run in sequence. The strict selection results and lax seletion results are derived from separate runs.

all_fgt_guides[, guides_per_tx := .N, by=transcript]
txs_per_gene = unique(all_fgt_guides[,.(transcript, guides_per_tx)])
all_gt = unique(coding_exons[, gene, by=transcript])

guides_per_tx = merge(all_gt, txs_per_gene, by=c("transcript"), all.x=TRUE)
guides_per_tx[is.na(guides_per_tx), guides_per_tx := 0]

#strict_selection_results = guides_per_tx
#strict_selection_results[, selection := "zero hg19 match at 0,1,2"]

#lax_selection_results = guides_per_tx
#lax_selection_results[, selection := "zero hg19 match at 0"]
 
#really_lax_selection_results = guides_per_tx
#really_lax_selection_results[, selection := "no hg19 restriction"]
 
supremely_lax_selection_results = guides_per_tx
supremely_lax_selection_results[, selection := "no hg19 restriction, zero mm10 match at 0,1"]

guides_improvement = data.table(rbind(strict_selection_results, lax_selection_results, really_lax_selection_results, supremely_lax_selection_results))
p1 = ggplot(guides_improvement, aes(x=guides_per_tx, fill=selection)) + geom_bar(position="dodge") + xlab("# guides per gene") + ylab("Number of genes")  +  ggtitle("Suitable guides per gene: 1371 genes total (all exons considered)")


