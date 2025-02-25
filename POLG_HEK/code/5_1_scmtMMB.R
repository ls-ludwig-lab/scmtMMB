setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
source("../../global_func/variant_calling.R")
source("../../global_func/quantifyMMB.R")
library(Matrix)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(data.table)
'%ni%' <- Negate('%in%')
color_vec <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")

# Load the seurat object
mtVar.bulk.df <- read.csv("../output/2_mtVar.bulk.meta.csv", row.names = "X")
CTRL_dnV <- mtVar.bulk.df %>% filter(occurrence == "de-novo" & condition == "CTRL") %>% pull(variant)
KI36_dnV <- mtVar.bulk.df %>% filter(occurrence != "Parental" & condition == "KI36") %>% pull(variant)
KIA2_dnV <- mtVar.bulk.df %>% filter(occurrence != "Parental" & condition == "KIA2") %>% pull(variant)

mmat.list <- list(readRDS("../POLG_HEK_large_data_files/output/2_mmtx.CTRL.rds")[CTRL_dnV,],
                  readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KI36.rds")[KI36_dnV,],
                  readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KIA2.rds")[KIA2_dnV,])
cov.list <- list(readRDS("../POLG_HEK_large_data_files/output/2_mtDNA_pos_coverage.CTRL.rds"),
                 readRDS("../POLG_HEK_large_data_files/output/2_mtDNA_pos_coverage.KI36.rds"),
                 readRDS("../POLG_HEK_large_data_files/output/2_mtDNA_pos_coverage.KIA2.rds"))

HEK_metadata <- read.csv("../output/1_HEK_metadata.csv",header = TRUE, row.names = 1)

scmtMMB_HEK <- do.call(rbind,lapply(1:3,function(idx){
  rbind(
    quantify_MMB(mmat.list[[idx]], cov.list[[idx]],names(color_vec)[idx]),
    quantify_MMB_psedobulk(mmat.list[[idx]], cov.list[[idx]],names(color_vec)[idx]))
}))

md <- HEK_metadata %>% dplyr::select(atacUMAP_1, atacUMAP_2,seurat_clusters, mtDNA_depth) %>% tibble::rownames_to_column("barcode")
scmtMMB_HEK <- scmtMMB_HEK %>% left_join(., md, by = "barcode")

write.csv(scmtMMB_HEK, "../output/5_scmtMMB_HEK.csv", quote = FALSE, row.names = TRUE)
