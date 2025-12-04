setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/deep_seq/")
source("../../global_func/variant_calling.R")
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(Matrix)
library(Signac)
library(Seurat)

mmat <- readRDS("output/2_mmtx.KI36.rds")
KI36 <- readRDS("output/1_POLG_HEK_preprocesssed_seurat.rds")[,colnames(mmat)]
KI36 <- AddMetaData(object = KI36,metadata = as.data.frame(KI36@reductions$atac.umap@cell.embeddings))
KI36 <- subset(KI36, atacUMAP_1 < 2 & atacUMAP_2 > 0)


var <-  as.data.frame(rowData(mmat)) %>% 
  subset(n_cells_conf_detected >= 50 &
                                                  n_cells_over_20 >50 &                          
                        strand_correlation > 0.65 &
                        vmr >= 0.01 &
                        variant %ni% c("301A>C", "302A>C", "310T>C", "316G>C")) %>% pull(variant)
mmat <- mmat[var,colnames(KI36)]


KI36[["allele"]] <- CreateAssayObject(assays(mmat)$allele_frequency)

DefaultAssay(KI36) <- "allele"
KI36 <- FindClonotypes(KI36, resolution = 0.02)
table(Idents(KI36))
DoHeatmap(KI36, features = VariableFeatures(KI36), slot = "data", disp.max = 0.1) +
  scale_fill_viridis_c()

DimPlot(KI36) + theme_void() + NoLegend()

variants_KI36 <- c("X15378T_C", "X13414G_A", "X10197G_A", "X8313G_A")
variant_names_KI36 <- c("m.15710C>T", "m.13414G>A", "m.10197G>A", "m.8313G>A")
VAF.long <- as.data.frame(KI36[["allele"]]$data) %>%
  rownames_to_column("variant") %>%
  separate(variant, into = c("Position", "Change"), sep = "(?<=\\d)(?=\\D)") %>%
  mutate(Position = as.numeric(Position)) %>%
  pivot_longer(cols = -c("Position", "Change"), names_to = "ID", values_to = "VAF")

VAF.long <- inner_join(VAF.long, md[,c("X","seurat_clusters", "condition", "atacUMAP_1", "atacUMAP_2")], by = c("ID" = "X"))
VAF.long <- VAF.long %>% mutate(present = VAF>0)

plot_mtSNV <- function(pos, var){
  
  p1 <- VAF.long %>% filter(Position == pos & Change == var) %>% arrange(VAF) %>% 
    ggplot(aes(atacUMAP_1, atacUMAP_2, color=VAF)) + geom_point() + scale_color_viridis_c(limits = c(0, 1), breaks = seq(0,1,0.2)) +
    ggtitle(paste0(pos, var)) + theme_void()
  
  #p2 <- VAF.long %>% filter(Position == pos & Change == var) %>% arrange(present) %>% 
  #  ggplot(aes(atacUMAP_1, atacUMAP_2, color=present)) + geom_point() + scale_color_manual(values = c("grey", "red")) + 
  #  theme_void() + ggtitle(paste0(pos, var))
  
  #p <- cowplot::plot_grid(p1, p2,ncol = 1)
  ggsave(plot = p1, paste0("../plot/Fig4_9_subcluster_clonal_var",pos, gsub(">", "", var),".pdf"), width = 2.5, height=2)
  p1
}

plot_mtSNV(7667, "C>T") # MT-CO2
plot_mtSNV(7927, "C>T") # MT-CO2
plot_mtSNV(15596, "G>A") # MT-CYB
plot_mtSNV(15378, "T>C") # MT-CYB
plot_mtSNV(10680, "G>A") # MT-ND4L
plot_mtSNV(15620, "C>T") # MT-ND4L
plot_mtSNV(6587, "C>T") # MT-CO2
plot_mtSNV(5660, "G>A") # MT-ND4L
