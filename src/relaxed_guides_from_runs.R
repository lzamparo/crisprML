library(data.table)

### Try to recover guides for fewer than covered and zero-coverage exons

### Load the guides from the main run and the needy exon run
setwd("~/projects/crisprML/data/")

### Try to rescue these genes by relaxing standards in human
rescue_txs = readRDS(file = "guides/needy_exon_run/rescue_txs.rds")

guide_features_mm10 = data.table(fread("guides/needy_exon_run/raw_features_computed_Cas9_sequences_vs_mm10.csv"))
guide_features_mm10 = guide_features_mm10[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_mm10 = unique(guide_features_mm10)

guide_features_hg19 = data.table(fread("guides/needy_exon_run/raw_features_computed_Cas9_sequences_vs_hg19.csv"))
guide_features_hg19 = guide_features_hg19[,.(sequence, Specificity_Score,Occurrences_at_Hamming_0,Occurrences_at_Hamming_1,Occurrences_at_Hamming_2,Occurrences_at_Hamming_3)]
guide_features_hg19 = unique(guide_features_hg19)

guides_to_targets = data.table(fread("guides/needy_exon_run/guide_to_target_region.csv", header=FALSE))
colnames(guides_to_targets) = c("guide", "chrom", "start", "end","strand")

coding_exons = data.table(fread("coding_exons/mm10_needy_exons.bed", sep="\t", header=TRUE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript", "transcript_start", "transcript_end", "gene")
coding_exons = unique(coding_exons)

### Join the guides based on sequence
setkey(guide_features_mm10, sequence)
setkey(guide_features_hg19, sequence)
feature_tables = merge(x=guide_features_hg19, y=guide_features_mm10, by.x="sequence", by.y="sequence", suffixes=c(".hg", ".mm"))
relaxed_hg_guides = feature_tables[Occurrences_at_Hamming_0.hg == 0 & Occurrences_at_Hamming_1.hg == 0 & Occurrences_at_Hamming_1.mm == 0 & Occurrences_at_Hamming_2.mm == 0,]
rm(guide_features_mm10, guide_features_hg19)

setkey(guides_to_targets, chrom, start, end)
setkey(coding_exons, chrom, start, end)
guides_to_exons = foverlaps(guides_to_targets, coding_exons, type="within", nomatch=0L)

### Start by assigning only those guides that are part of the filtered guides set
setkey(relaxed_hg_guides, sequence)
setkey(guides_to_exons, guide)
filtered_guides_and_targets = merge(relaxed_hg_guides, guides_to_exons, by.x=c("sequence"), by.y=c("guide"))
setnames(filtered_guides_and_targets, c("i.start", "i.end"), c("guide_start", "guide_end"))


### select those guides which target rescue_txs
filtered_guides_and_targets = filtered_guides_and_targets[transcript %in% rescue_txs[,transcript],]
filtered_guides_and_targets[, num_guides := uniqueN(sequence), by = .(gene)]
filtered_guides_and_targets[, exon_count := uniqueN(exon), by = gene]
filtered_guides_and_targets[, hamming_3_sum := Occurrences_at_Hamming_3.mm + Occurrences_at_Hamming_3.hg]

### Break down filtered_guides_and_targets by gene with 4,3,2,1 exons
four_exons_fgt = filtered_guides_and_targets[exon_count >= 4,]
three_exons_fgt = filtered_guides_and_targets[exon_count == 3,]
two_exons_fgt = filtered_guides_and_targets[exon_count == 2,]
one_exons_fgt = filtered_guides_and_targets[exon_count == 1,]

### Four: choose one guide / exon, each with minimum hamming_3_sum, first four exons by position
four_exon_fgt_guides = four_exons_fgt[, .SD[which.min(hamming_3_sum)], by = .(exon)]
four_exon_fgt_guides = four_exon_fgt_guides[,.(sequence, exon, transcript, gene, chrom, start, end, hamming_3_sum, guide_start, guide_end)]
four_exon_fgt_guides = four_exon_fgt_guides[,.SD[head(order(start),4)], by = .(transcript)]

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

### Combine all guides together, write out guides
all_relaxed_guides = data.table(rbind(four_exon_fgt_guides,three_exon_fgt_guides,two_exon_fgt_guides,one_exon_fgt_guides))
write.csv(all_relaxed_guides,file = "../results/library/2531_relaxed_guides.csv", row.names=FALSE)



