library(data.table)
library(ggplot2)
library(biomaRt)

### extract the list of transcripts targeted by Doench in the Brie library, 
### keep only those exons belonging to that set, write out the
### set of first two or first four coding exons of each

setwd("/Users/zamparol/projects/crisprML/data")

### Get our table of first four coding exons for all txs, and the Brie list of txs
coding_exons = data.table(fread("coding_exons/mm10_first_four_coding_exons_ensembl.bed", sep="\t", header=TRUE))
colnames(coding_exons) = c("chrom","start","end","strand","exon","transcript","gene", "transcript_start", "transcript_end")

doench_guides = data.table(fread("coding_exons/broadgpp-brie-library-contents.txt", sep="\t", header=TRUE))
### guard against different defaults in fread introducing spaces (or not replacing them with '.')
if ("Target Transcript" %in% colnames(doench_guides)){
  setnames(doench_guides, "Target Transcript", "Target.Transcript")
}
doench_guides[, transcript := tstrsplit(Target.Transcript, ".", fixed=TRUE)[[1]]]

### get the ENSEMBL transcript codes for the Brie RefSeq txs
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
txs = as.character(unlist(doench_guides[!duplicated(transcript),transcript]))
res = getBM(attributes=c("refseq_mrna", "ensembl_transcript_id"),
            filters = c("refseq_mrna"), 
            values = txs, 
            mart = mouse)
res_dt = data.table(res)


### keep only the exons int the res_dt txs set
setkey(res_dt, ensembl_transcript_id)
setkey(coding_exons, transcript)
filtered_txs = coding_exons[transcript %in% unique(res_dt[, ensembl_transcript_id]),]

### keep only the relevant stuff
write.table(filtered_txs, file = "coding_exons/mm10_first_four_filtered_brie_exons.bed", sep="\t", quote=FALSE, row.names=FALSE)
