library(data.table)
library(ggplot2)
library(stringr)

### Load the guides from the main run and the needy exon run
setwd("~/projects/crisprML/results/library/")
main_run_guides = data.table(read.csv(file="main_run_four_guides_per_gene_final.csv"))
combined_rescue = data.table(read.csv(file="rescue_guides_for_low_genes.csv"))

#needy_exon_run_guides = data.table(read.csv(file="supplementary_run_local.csv"))
#missing_exon_run_guides = data.table(read.csv(file="missing_exon_run.csv"))
#all_fgt_guides = data.table(rbind(main_run_guides, needy_exon_run_guides[,.(sequence, exon, transcript, gene, chrom, start, end, hamming_3_sum, guide_start, guide_end)], missing_exon_run_guides[,.(sequence, exon, transcript, gene, chrom, start, end, hamming_3_sum, guide_start, guide_end)]))

all_fgt_guides = unique(data.table(rbind(main_run_guides, combined_rescue)))
setkey(all_fgt_guides, transcript, exon)

### Load the map of all exons to be covered
setwd("~/projects/crisprML/data/")
coding_exons = data.table(fread("coding_exons/mm10_first_four_filtered_brie_exons.bed", sep="\t", header=TRUE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript","gene", "transcript_start", "transcript_end")
coding_exons = unique(coding_exons)

### Which genes do not yet meet the four guides / exon threshold?
### How many still do not have four guides / gene?
all_fgt_guides[, exon_count := uniqueN(exon), by = gene]
four_exons_fgt = all_fgt_guides[exon_count >= 4,]
three_exons_fgt = all_fgt_guides[exon_count == 3,]
two_exons_fgt = all_fgt_guides[exon_count == 2,]
one_exons_fgt = all_fgt_guides[exon_count == 1,]

### Four: choose one guide / exon, each with minimum hamming_3_sum, first four exons by position
four_exon_fgt_guides = four_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
#four_exon_fgt_guides = four_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]
four_exon_fgt_guides = four_exon_fgt_guides[,.SD[head(order(start),4)], by = .(transcript)]

### Three: choose one guide / exon, and then the minimum hamming_3_sum of remaining guides
three_exon_fgt_guides = three_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
chosen_guides = three_exon_fgt_guides[,sequence]
additionals = three_exons_fgt[!(sequence %in% chosen_guides), .SD[head(order(hamming_3_sum),1)], by = .(transcript)]
three_exon_fgt_guides = data.table(rbind(three_exon_fgt_guides,additionals))
#three_exon_fgt_guides = three_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]

### Two: choose one guide / exon, and then the smallest two hamming_3_sum of remaining guides
two_exon_fgt_guides = two_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
chosen_guides = two_exon_fgt_guides[,sequence]
additionals = two_exons_fgt[!(sequence %in% chosen_guides), .SD[head(order(hamming_3_sum),2)], by = .(transcript)]
two_exon_fgt_guides = data.table(rbind(two_exon_fgt_guides,additionals))
#two_exon_fgt_guides = two_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]


### One: choose the four guides with minimal hamming_3_sum
one_exon_fgt_guides = one_exons_fgt[, .SD[head(order(hamming_3_sum),4)], by = .(exon)]
#one_exon_fgt_guides = one_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, Occurrences_at_Hamming_2.mm, Occurrences_at_Hamming_2.hg, Occurrences_at_Hamming_1.hg, Occurrences_at_Hamming_0.hg, hamming_3_sum, guide_start, guide_end)]

### combine the data from all runs 
all_fgt_guides = data.table(rbind(four_exon_fgt_guides,three_exon_fgt_guides,two_exon_fgt_guides,one_exon_fgt_guides))
all_fgt_guides[, num_guides := uniqueN(sequence), by=transcript]

### test the set to make sure, one last time, there are no damn spurious matches for the restriction enzymes

# BstXI: 73 guides need fixing
# most seem to be barcode problem (see below)
# others seem to be misses on existing filters (maybe these are controls?  They are controls.)

# Barcode problem: C CACGTTCCTGG  don't want any BCs matching ^CA[ACGT]{6}TGG

BamHI_filter_emergent = "^[G]?ATCC"
BlpI_filter_emergent = "^CT[ACGT]{1}AGC"
BstXI_filter_emergent = "CCA[ACGT]{6}TG$"

BstXI_filter_exact = "CCA[ACGT]{6}TGG"
BamHI_filter_exact = "GGATCC"
BlpI_filter_exact = "GCT[ACGT]{1}AGC"

BstXI_filter_barcode = "^CA[ACGT]{6}TGG"

all_fgt_guides[str_count(sequence, BstXI_filter_exact) > 0,.N]
all_fgt_guides[str_count(sequence, BamHI_filter_exact) > 0,.N]
all_fgt_guides[str_count(sequence, BlpI_filter_exact) > 0,.N]
all_fgt_guides[str_count(sequence, BamHI_filter_emergent) > 0,.N]
all_fgt_guides[str_count(sequence, BstXI_filter_emergent) > 0,.N]
all_fgt_guides[str_count(sequence, BlpI_filter_emergent) > 0,.N]

### Write out the guides for genes with enough coverage
write.csv(all_fgt_guides, file="../results/library/gene_library.csv", row.names=FALSE)

