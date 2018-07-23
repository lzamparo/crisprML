library(data.table)
library(readxl)

### Read and assemble the consructs from gene_library

# load the data
curr_dir = getwd()
setwd("~/projects/crisprML/results/library/")

guides_dt = data.table(fread("gene_library.csv"))
provided_xls = as.data.table(read_xlsx("Machine Learning Library-Master File.xlsx"))
barcodes = data.table(read.csv("../new_barcodes.txt", header = FALSE))
colnames(barcodes) = c("barcode")

# regular expressions to find and eliminate spurious matches to the restriction enzymes
BamHI_filter_emergent = "^[G]?ATCC"
BlpI_filter_emergent = "^CT[ACGT]{1}AGC"
BstXI_filter_emergent = "CCA[ACGT]{6}TG$"
BstXI_filter_barcode = "^CA[ACGT]{6}TGG"

BstXI_filter_exact = "CCA[ACGT]{6}TGG"
BamHI_filter_exact = "GGATCC"
BlpI_filter_exact = "GCT[ACGT]{1}AGC"

# test to make sure the constructs are merged in order of L -> R columns
provided_xls[, test_order := paste0(Forward,BstxI,Guides,Scaffold,Barcode,Buffer,Target,PAM,BamHI,Reverse)]
stopifnot(provided_xls[!is.na(Target), .N] == provided_xls[!is.na(Target) & (test_order == Merged), .N])
provided_xls[,test_order := NA]

# drop rows that have no guide
provided_xls = provided_xls[!is.na(Target),]

# replace guides for non-controls with guides
guides = guides_dt[,sequence]
guide_names = guides_dt[, gene]
provided_xls[is.na(`Gene ID`),Target := guides]
provided_xls[is.na(`Gene ID`), `Gene ID` := guide_names]

# find & eliminate controls that have spurious matches for the restriction enzymes in the guides
provided_xls = provided_xls[str_count(Target, BstXI_filter_exact) == 0,]
provided_xls = provided_xls[str_count(Target, BamHI_filter_exact) == 0,]
provided_xls = provided_xls[str_count(Target, BlpI_filter_exact) == 0,]
provided_xls = provided_xls[str_count(Target, BamHI_filter_emergent) == 0,]
provided_xls = provided_xls[str_count(Target, BstXI_filter_emergent) == 0,]
provided_xls = provided_xls[str_count(Target, BlpI_filter_emergent) == 0,]

# find & replace barcodes for guides with partial match for BstXI's recognition site
used_barcodes = provided_xls[,Barcode]
remaining_barcodes = barcodes[!(barcode %in% used_barcodes),]
remaining_barcodes = remaining_barcodes[(str_count(barcode, BstXI_filter_barcode) == 0) & (str_count(barcode, BstXI_filter_emergent) == 0),]
remaining_barcodes = remaining_barcodes[(str_count(barcode,BamHI_filter_exact) == 0) & (str_count(barcode,BamHI_filter_emergent) == 0),]
remaining_barcodes = remaining_barcodes[(str_count(barcode,BlpI_filter_exact) == 0) & (str_count(barcode,BlpI_filter_emergent) == 0),]

replace_these_codes = provided_xls[str_count(Barcode, BstXI_filter_barcode) > 0 | str_count(Barcode, BstXI_filter_emergent) > 0 | str_count(Barcode,BamHI_filter_exact) > 0 | str_count(Barcode,BamHI_filter_emergent) > 0 | str_count(Barcode,BlpI_filter_exact) > 0 | str_count(Barcode,BlpI_filter_emergent) > 0,Barcode]
replacement_codes = remaining_barcodes[sample(.N,length(replace_these_codes))]
provided_xls[Barcode %in% replace_these_codes, Barcode := replacement_codes]

# replace the constructs with updated barcodes
provided_xls[, Merged := paste0(Forward,BstxI,Guides,Scaffold,Barcode,Buffer,Target,PAM,BamHI,Reverse)]

# write out construct library
setnames(provided_xls, "Merged", "Constructs")
write.csv(provided_xls, "construct_library.csv", row.names=FALSE,quote=FALSE)
