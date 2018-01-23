library(data.table)
library(stringr)

### Find those constructs with spurious matches of the restriction enzymes BstxI, BamHI, in the constructs

# load the constructs
setwd('~/projects/crisprML/results/library/')
construct_library = data.table(fread("library_sequences_fixed.csv"))

# what are the motifs for the restriction enzymes?
BstXI_motif = "CCANNNNNNTGG"
BamHI_motif = "GGATCC"
BlpI_motif = "GCTNAGC"

BstXI_match = "CCA[ACGT]{6}TGG"
BamHI_match = BamHI_motif
BlpI_match = "GCT[ACGT]{1}AGC"

# count occurrences in the barcodes
construct_library[,BstXI_barcode_count := str_count(Barcode,BstXI_match)]
construct_library[,BamHI_barcode_count := str_count(Barcode,BamHI_match)]
construct_library[,BlpI_barcode_count := str_count(Barcode,BlpI_match)]

# count occurrences in the guide (Guides in the library)
construct_library[,BstXI_guides_count := str_count(Guides,BstXI_match)]
construct_library[,BamHI_guides_count := str_count(Guides,BamHI_match)]
construct_library[,BlpI_guides_count := str_count(Guides,BlpI_match)]

# count occurrences outside between the scaffold, barcode and buffer
construct_library[,BstXI_construct_count := str_count(str_to_upper(Guide),BstXI_match) - str_count(Guides,BstXI_match) - str_count(Barcode,BstXI_match)]
construct_library[,BamHI_construct_count := str_count(str_to_upper(Guide),BamHI_match) - str_count(Guides,BamHI_match) - str_count(Barcode,BamHI_match)]
construct_library[,BlpI_construct_count := str_count(str_to_upper(Guide),BlpI_match) - str_count(Guides,BlpI_match) - str_count(Barcode,BlpI_match)]

# find those guides or barcodes which need replacing: guide or barcode count > 0
guide_or_barcode_count_violations = unique(construct_library[BstXI_barcode_count > 0 | BamHI_barcode_count > 0 | BlpI_barcode_count > 0 | BstXI_guides_count > 1 | BamHI_guides_count > 1 | BlpI_guides_count > 1, .(Guide, Guides, Barcode)])

# find those guides or barcodes which need replacing: construct counts > 1
construct_counts_violations = unique(construct_library[BstXI_construct_count > 1 | BamHI_construct_count > 1 | BlpI_construct_count > 1, .(Guide, Guides, Barcode)])

# merge the two sets
constructs_needing_work = unique(data.table(rbind(guide_or_barcode_count_violations, construct_counts_violations)))

# find where the extra matches occur



# for those constructs for matches in the Txsfind which Txs 


