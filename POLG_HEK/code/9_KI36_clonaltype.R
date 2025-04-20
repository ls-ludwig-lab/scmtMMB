setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code")
source("../../global_func/variant_calling.R")
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(Matrix)
library(Signac)
library(Seurat)

mmat <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KI36.rds")
KI36 <- readRDS("../POLG_HEK_large_data_files/output/1_POLG_HEK_preprocesssed_seurat.rds")[,colnames(mmat)]
KI36 <- AddMetaData(object = KI36,metadata = as.data.frame(KI36@reductions$atac.umap@cell.embeddings))
KI36 <- subset(KI36, atacUMAP_1 <=4 & atacUMAP_2 >=2)


var <-  as.data.frame(rowData(mmat)) %>% 
  subset(n_cells_conf_detected >= 50 &
                                                  n_cells_over_20 >50 &                          
                        strand_correlation > 0.65 &
                        vmr >= 0.01 &
                        variant %ni% c("301A>C", "302A>C", "310T>C", "316G>C")) %>% pull(variant)
mmat <- mmat[var,colnames(KI36)]


KI36[["allele"]] <- CreateAssayObject(assays(mmat)$allele_frequency)

DefaultAssay(KI36) <- "allele"
KI36 <- FindClonotypes(KI36, resolution = 0.3)
table(Idents(KI36))
DoHeatmap(KI36, features = VariableFeatures(KI36), slot = "data", disp.max = 0.1) +
  scale_fill_viridis_c()

DimPlot(KI36)

px2 <- DoHeatmap(KI36, slot = "data", disp.max = 1, 
                 features = as.data.frame(rowData(mmat)) %>% 
                   arrange(desc(mean), desc(n_cells_conf_detected)) %>%
                   slice_head(n = 10) %>%  pull(variant)) + scale_fill_viridis_c() + NoLegend()
