library(data.table)
library(stringr)

### Find those constructs with spurious matches of the restriction enzymes BstxI, BamHI, in the constructs

# load the constructs
setwd('~/projects/crisprML/results/library/')
construct_library = data.table(fread("construct_library.csv"))

# what are the motifs for the restriction enzymes?
BamHI_filter_emergent = "^[G]?ATCC"
BlpI_filter_emergent = "^CT[ACGT]{1}AGC"
BstXI_filter_emergent = "CCA[ACGT]{6}TG$"
BstXI_filter_barcode = "^CA[ACGT]{6}TGG"

BstXI_filter_exact = "CCA[ACGT]{6}TGG"
BamHI_filter_exact = "GGATCC"
BlpI_filter_exact = "GCT[ACGT]{1}AGC"

BstXI_filter_barcode = "^CA[ACGT]{6}TGG"

# count occurrences in the barcodes
construct_library[,BstXI_barcode_count := str_count(Barcode,BstXI_filter_exact) + str_count(Barcode,BstXI_filter_emergent) + str_count(Barcode,BstXI_filter_barcode)]
construct_library[,BamHI_barcode_count := str_count(Barcode,BamHI_filter_exact) + str_count(Barcode,BamHI_filter_emergent)]
construct_library[,BlpI_barcode_count := str_count(Barcode,BlpI_filter_exact) + str_count(Barcode,BlpI_filter_emergent)]

# count occurrences in the guide (Guides in the library)
construct_library[,BstXI_guides_count := str_count(Guides,BstXI_filter_exact) + str_count(Guides,BstXI_filter_emergent)]
construct_library[,BamHI_guides_count := str_count(Guides,BamHI_filter_exact) + str_count(Guides,BamHI_filter_emergent)]
construct_library[,BlpI_guides_count := str_count(Guides,BlpI_filter_exact) + str_count(Guides,BlpI_filter_emergent)]

# validate barcode constraints
stopifnot(construct_library[BstXI_barcode_count > 0,.N] == 0)
stopifnot(construct_library[BamHI_barcode_count > 0,.N] == 0)
stopifnot(construct_library[BlpI_barcode_count > 0,.N] == 0)

# validate guides constraints
stopifnot(construct_library[BstXI_guides_count > 0,.N] == 0)
stopifnot(construct_library[BamHI_guides_count > 0,.N] == 0)
stopifnot(construct_library[BlpI_guides_count > 0,.N] == 0)

# count occurrences outside between the guide, scaffold, barcode and buffer
construct_library[,BstXI_construct_count := str_count(str_to_upper(Constructs), BstXI_filter_exact)]
construct_library[,BamHI_construct_count := str_count(str_to_upper(Constructs), BamHI_filter_exact)]
construct_library[,BlpI_construct_count := str_count(str_to_upper(Constructs), BlpI_filter_exact)]

stopifnot(construct_library[BstXI_construct_count > 1, .N] == 0)
stopifnot(construct_library[BamHI_construct_count > 1, .N] == 0)
stopifnot(construct_library[BlpI_construct_count > 1, .N] == 0)

