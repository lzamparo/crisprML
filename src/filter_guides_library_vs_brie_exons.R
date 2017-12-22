### Characterize a given set of exons wrt how many guides they have per gene, of what quality, etc.
library(data.table)
library(ggplot2)

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
filtered_guides = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_1.hg == 0 & Occurrences_at_Hamming_2.hg == 0 & Occurrences_at_Hamming_0.mm == 1 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]

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
filtered_guides_and_targets[, num_guides := .N, by = .(gene)]
genes_and_guides = unique(filtered_guides_and_targets[,.(gene, num_guides)])
covered_genes = genes_and_guides[num_guides >= 4,gene]
filtered_guides_and_targets[, exon_count := uniqueN(exon), by = gene]

### What about just the guide which targets each exon, but with minimal Hamming_3.mm + Hamming_3.hg?
filtered_guides_and_targets[, hamming_3_sum := Occurrences_at_Hamming_3.mm + Occurrences_at_Hamming_3.hg]

### Break down filtered_guides_and_targets by gene with 4,3,2,1 exons
four_exons_fgt = filtered_guides_and_targets[gene %in% covered_genes & exon_count == 4,]
three_exons_fgt = filtered_guides_and_targets[gene %in% covered_genes & exon_count == 3,]
two_exons_fgt = filtered_guides_and_targets[gene %in% covered_genes & exon_count == 2,]
one_exons_fgt = filtered_guides_and_targets[gene %in% covered_genes & exon_count == 1,]


### Four: choose one guide / exon, each with minimum hamming_3_sum
four_exon_fgt_guides = four_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]

### Three: choose one guide / exon, and then the minimum hamming_3_sum of remaining guides
three_exon_fgt_guides = three_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
chosen_guides = three_exon_fgt_guides[,sequence]
additionals = three_exons_fgt[!(sequence %in% chosen_guides), .SD[head(order(hamming_3_sum),1)], by = .(transcript)]
three_exon_fgt_guides = 

### Two: choose one guide / exon, and then the smallest two hamming_3_sum of remaining guides
two_exon_fgt_guides = two_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]

### One: choose the four guides with minimal hamming_3_sum
one_exon_fgt_guides = one_exons_fgt[, .SD[head(order(hamming_3_sum),4)], by = .(exon)]
#d[,.SD[head(order(mpg))],by=cyl]


min_off_target_guides_per_exon = filtered_guides_and_targets[, .SD[which.min(hamming_3_sum)], by = .(exon)]
ggplot(min_off_target_guides_per_exon, aes(x=hamming_3_sum)) + geom_histogram(bins=100) + xlim(0,100) + xlab("Number of matches with up to 3 mismatches in hg19 + mm10") + ylab("# guides") + ggtitle("Histogram of Brie exome guide specificity")

### Select the min off-target guides that target exons in covered genes
good_guides = min_off_target_guides_per_exon[gene %in% covered_genes,.(sequence,chrom,start,end,strand,exon,transcript,gene,hamming_3_sum)]


# number of suitable guides per exon
g_per_e = unique(filtered_guides_and_targets[order(-num_guides),.(gene, transcript, num_guides),by=exon])

# (1) those exons with suitable guides
exons_w_suitable_guide = g_per_e[,exon]

# those exons with at least one guide, suitable or not
exons_w_guide = guides_to_exons[,.N,by=exon][order(-N),exon]

# (2) those exons without a suitable guide
exons_wo_suitable_guide = setdiff(exons_w_guide, exons_w_suitable_guide)

# (3) those exons with no guide at all
exons_w_no_guides = coding_exons[!(exon %in% exons_w_guide),]

# (1) + (2) + (3) = total exon count 
dim(g_per_e)[1] + length(exons_wo_suitable_guide) + dim(exons_w_no_guides)[1]

# again, to make sure
all_exons = c(exons_w_no_guides[,exon], exons_wo_suitable_guide, exons_w_suitable_guide)
length(all_exons)
# account for multiplicity by gene synonyms, and one shared exons in Foxa1
length(unique(all_exons))

