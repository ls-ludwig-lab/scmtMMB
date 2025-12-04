setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/deep_seq/")

# Load the seurat object
HEK <- readRDS("../POLG_HEK_large_data_files/output/1_POLG_HEK_preprocesssed_seurat.rds")

# extract cell barcodes for each line
CTRL_CB <- WhichCells(HEK, idents = "CTRL")
KI36_CB <- WhichCells(HEK, idents = "KI36")
KIA2_CB <- WhichCells(HEK, idents = "KIA2")


writeLines(paste0("CB:Z:", CTRL_CB), "CTRL_CB.txt")
writeLines(paste0("CB:Z:", KI36_CB), "KI36_CB.txt")
writeLines(paste0("CB:Z:", KIA2_CB), "KIA2_CB.txt")


writeLines(CTRL_CB, "CTRL_CB_mgatk.tsv")
writeLines(KI36_CB, "KI36_CB_mgatk.tsv")
writeLines(KIA2_CB, "KIA2_CB_mgatk.tsv")
