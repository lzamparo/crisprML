library(data.table)

### Get exon, transcript, gene name info from guides from GUIDES which have no coverage in the present library
setwd("/Users/zamparol/projects/crisprML/data")

# load the coding exons
coding_exons = data.table(fread("coding_exons/mm10_first_four_filtered_brie_exons.bed", sep="\t", header=TRUE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript","gene", "transcript_start", "transcript_end")

setwd("/Users/zamparol/projects/crisprML/results/library")
# load the existing library
existing_lib = data.table(fread("permissive_18995_gene_library.csv"))

# select those txs which with low 
needy_txs = existing_lib[num_guides < 4,.(gene,transcript)]

# (or no) coverage
txs_list = unique(existing_lib[,transcript])
zero_coverage_txs = coding_exons[!(transcript %in% txs_list),.(gene,exon,transcript,transcript_start,transcript_end)]

# merge gene list from needy_txs, zero_coverage_txs, query GUIDES for the guides to cover these (offline)
my_genes = unique(data.table(rbind(needy_txs,zero_coverage_txs[,.(gene,transcript)])))
write.table(my_genes[,gene],"/Users/zamparol/projects/crisprML/data/guides/guides_from_GUIDES/genes.txt",sep="\n",col.names=FALSE,row.names=FALSE,quote=FALSE)
my_results = data.table(fread("/Users/zamparol/projects/crisprML/data/guides/guides_from_GUIDES/guides.csv"))

# clean up results, merge with coding_exons by foverlaps to get description of which exons guides are targeting
my_results[,guide_start := pam_guide_end - 20]
my_results[,chrom := paste0("chr",boguschrom)]
setkey(my_results,chrom, guide_start,pam_guide_end)
setkey(coding_exons,chrom,start,end)

overlaps_w_needy_txs = foverlaps(my_results,coding_exons,by.x=c("chrom","guide_start","pam_guide_end"),by.y=c("chrom","start","end"),type="within",nomatch=NA)
needy_txs_guides = overlaps_w_needy_txs[!is.na(start),.(sequence,exon,transcript,gene,chrom,start,end,transcript_start,transcript_end,pam_guide_end)]
needy_txs_guides[, num_guides := .N, by=transcript]
needy_txs_guides[, hamming_3_sum := NA]
needy_txs_guides[, guide_start := pam_guide_end - 20]
setnames(needy_txs_guides, "pam_guide_end", "guide_end")

# merge with existing guides to get final library
supplemental_guides = needy_txs_guides[,.(sequence,exon,transcript,gene,chrom,start,end,hamming_3_sum,guide_start,guide_end,num_guides)]
supplemental_guides[,gs_guide := FALSE]
existing_lib[,gs_guide := TRUE]

merged_library = data.table(rbind(existing_lib,supplemental_guides))
write.csv(merged_library, file="merged_permissive_library_w_GUIDES.csv", row.names=FALSE)
