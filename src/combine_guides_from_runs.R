library(data.table)
library(ggplot2)

### Load the guides from the main run and the needy exon run
setwd("~/projects/crisprML/results/library/")
main_run_guides = data.table(read.csv(file="main_run_four_guides_per_gene.csv"))
needy_exon_run_guides = data.table(read.csv(file="supplementary_run.csv"))
all_fgt_guides = data.table(rbind(main_run_guides, needy_exon_run_guides[,.(sequence, exon, transcript, gene, chrom, start, end, hamming_3_sum, guide_start, guide_end)]))
setkey(all_fgt_guides, transcript, exon)

### Load the map of all exons to be covered
setwd("~/projects/crisprML/data/")
coding_exons = data.table(fread("coding_exons/mm10_first_four_filtered_brie_exons.bed", sep="\t", header=TRUE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript","gene", "transcript_start", "transcript_end")
coding_exons = unique(coding_exons)

### Which genes do not yet meet the four guides / exon threshold?
### How many still do not have four guides / gene?
guides_per_tx = unique(all_fgt_guides[, num_guides := .N, by = transcript])
guides_per_tx = unique(guides_per_tx[,.(transcript, num_guides)])

### Are there any genes yet to be covered?
setkey(coding_exons, transcript)
setkey(guides_per_tx, transcript)
coverage = merge(guides_per_tx, coding_exons, by.x="transcript", by.y="transcript", all.y=TRUE)

### Fewer than 4 guides / tx
under_covered_txs = coverage[!is.na(num_guides) & num_guides < 4, ]
uc_txs = unique(under_covered_txs[,.(chrom,transcript_start,transcript_end,gene,num_guides),by=transcript])

### No guides / tx
no_coverage = coverage[is.na(num_guides),]
nc_txs = unique(no_coverage[,.(chrom,transcript_start,transcript_end,gene),by=transcript])
nc_txs[,num_guides := 0]

### Try to rescue these genes by relaxing standards in human
rescue_txs = unique(data.table(rbind(uc_txs,nc_txs)))
saveRDS(rescue_txs,file = "guides/needy_exon_run/rescue_txs.rds")

### Which genes do not yet meet the four guides / exon threshold?
### How many still do not have four guides / gene?
guides_per_tx = unique(all_fgt_guides[, num_guides := .N, by = transcript])
guides_per_tx = unique(guides_per_tx[,.(transcript, num_guides)])
all_fgt_guides[num_guides >= 4,]
write.csv(all_fgt_guides[num_guides >= 4,], file="../results/library/conservative_18100_gene_library.csv", row.names=FALSE)

### add relaxed guides?
all_relaxed_guides = data.table(read.csv(file = "../results/library/2531_relaxed_guides.csv"))
all_relaxed_guides[, num_guides := .N, by = transcript]
all_combined_guides = unique(data.table(rbind(all_fgt_guides,all_relaxed_guides)))
write.csv(all_combined_guides, file="../results/library/permissive_18995_gene_library.csv", row.names=FALSE)
