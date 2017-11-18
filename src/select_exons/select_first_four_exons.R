library(data.table)
library(biomaRt)

#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")

# set up connection to biomaRt
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
listFilters(mouse)

# get all protein coding exons, annotated with transcript id, start, end, etc.
res = getBM(attributes=c("ensembl_transcript_id","transcript_start","transcript_end","ensembl_exon_id","exon_chrom_start","exon_chrom_end","strand","chromosome_name","gene_biotype"),
            filters = c("chromosome_name","biotype"), 
            values=list(c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y"),"protein_coding"),
            mart=mouse)
res_dt = data.table(res)

# now get all gene names associated with each transcript
gene_names = getBM(attributes=c("ensembl_transcript_id","mgi_symbol","chromosome_name"),
                    filters = c("chromosome_name"), 
                    values=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y"),
                    mart=mouse)
gene_names_dt = data.table(gene_names)
gene_names_dt = gene_names_dt[,.(ensembl_transcript_id, mgi_symbol)]

# merge res_dt with gene_names
setkey(res_dt, ensembl_transcript_id)
setkey(gene_names_dt, ensembl_transcript_id)
coding_exons_gene_names_dt = merge(res_dt, gene_names_dt, by.x=c("ensembl_transcript_id"), by.y=c("ensembl_transcript_id"))

first_four_coding_exons_by_tx_dt = coding_exons_gene_names_dt[,  head(.SD[order(exon_chrom_start)],4), by=.(ensembl_transcript_id)]

#d <- data.table(mtcars, key="cyl")
#d[, head(.SD[order(disp)], 3), by=cyl]

grouped_filtered_guides_to_remaining_dt = filtered_guides_to_remaining_exons[, .SD[which.min(score_summary)], by=.(gene, transcript, exon_number)]


exons = data.table(read.delim("mm10_annotation_cds_only_coding_exons_sorted_unique_bygene_exon.bed", sep="\t", header = FALSE))
colnames(exons) = c("chrom", "start", "end", "gene", "transcript", "exon_num", "max_exons")

first_four_exons <- exons %>% group_by(gene) %>% top_n(-4, exon_num)
first_four = as.data.table(ungroup(first_four_exons))

ff_exons = as.data.table(ungroup(first_five_exons))
write.table(ff_exons, file="mm10_first_five_exons.bed", sep="\t", row.names=FALSE, quote=FALSE)
