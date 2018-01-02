library(ggplot2)
library(data.table)

### Some measure of validation: how often do we recover our guides when using the GUIDES tool?

setwd("/Users/zamparol/projects/crisprML/results/guides_tests")

# read in the list of genes, and the dt of guides
gene_list_dt = data.table(read.table(file = "gene_list.txt", sep="\n"))
colnames(gene_list_dt) = c("genes")
our_guides = data.table(fread("/Users/zamparol/projects/crisprML/results/library/permissive_18995_gene_library.csv"))

# sample without replacement the list of genes: 10 samples, 500 genes each

for (j in 1:10){
  dirname = paste("sample",j,sep="_")
  sample_w_regions_name = "sample_genes_w_regions.csv"
  sample_w_genes_name = "sample_genes.txt"
  if (!dir.exists(dirname)){
    dir.create(dirname)
  }
  gene_sample = gene_list_dt[sample(1:18995,500,replace=FALSE),]
  sample_w_regions = our_guides[gene %in% gene_sample$genes,.(sequence,gene,chrom,start,end)]
  write.csv(sample_w_regions, file=paste(dirname,sample_w_regions_name,sep="/"),row.names=FALSE)
  write.table(gene_sample, file=paste(dirname,sample_w_genes_name,sep="/"),row.names=FALSE,col.names=FALSE,sep="\n")
}


# (out of script): get the list of guides from GUIDES
### done
overlap_count = c()
match_count = c()

# Calculate interesection per sample of guides for genes
for (j in 1:10){
  guides_sample = data.table(fread(paste(paste0("sample_",j),"guides.csv",sep="/")))
  guides_chosen_by_us = data.table(fread(paste(paste0("sample_",j),"sample_genes_w_regions.csv",sep="/")))
  setnames(guides_sample, c("Chromosome","PAM position in chromosome","Targets last exon"),c("chrom","guide_end","targets_last_exon"))
  guides_sample[, guide_start := guide_end - 20]
  guides_sample[, adj_chrom := paste0("chr",chrom)]
  setkey(guides_sample, adj_chrom, guide_start, guide_end)
  setkey(guides_chosen_by_us, chrom, start, end)
  
  # how many guides for the same gene are shared?
  # how many guides for the same genes targeting the same exons are shared?
  joined_guides = foverlaps(guides_sample, guides_chosen_by_us, by.x=c("adj_chrom","guide_start", "guide_end"), by.y=c("chrom","start","end"), type="within", nomatch=NA, mult="first")
  joined_guides = joined_guides[!is.na(sequence),]
  num_overlaps = joined_guides[Sequence == sequence, .N]
  overlap_count = c(overlap_count, num_overlaps)
  match_count = c(match_count, joined_guides[,.N])
}

