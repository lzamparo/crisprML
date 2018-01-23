library(data.table)
library(ggplot2)
library(biomaRt)

### extract the list of transcripts targeted by Doench in the Brie library, 
### keep only those exons belonging to that set, write out the
### set of first two or first four coding exons of each

setwd("/Users/zamparol/projects/crisprML/data")

### get the list of Txs needing more guides
txs_needing_more_guides = data.table(fread("coding_exons/txs_needing_more_guides.csv", header=TRUE))

### get the ENSEMBL transcript codes for the txs with low (or no) coverage
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
txs = as.character(unlist(unique(txs_needing_more_guides[,transcript])))
res = getBM(attributes=c("ensembl_exon_id", "ensembl_transcript_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", "external_gene_name", "transcript_start", "transcript_end"),
            filters = c("ensembl_transcript_id"), 
            values = txs, 
            mart = mouse)
res_dt = data.table(res)

### keep only the exons in the res_dt txs set
### format output like: chrom	start	end	strand	exon	transcript	gene	transcript_start	transcript_end
setnames(res_dt, c("chromosome_name", "exon_chrom_start", "exon_chrom_end", "external_gene_name", "ensembl_transcript_id", "ensembl_exon_id"), c("chrom", "start","end","gene","transcript","exon"))
res_dt[,mychrom := paste0("chr",chrom)]
filtered_txs = res_dt[,.(mychrom, start, end, strand, exon, transcript, transcript_start, transcript_end, gene)]
setnames(filtered_txs, "mychrom", "chrom")

### keep only the relevant stuff
write.table(filtered_txs, file = "coding_exons/mm10_needy_exons.bed", sep="\t", quote=FALSE, row.names=FALSE)
