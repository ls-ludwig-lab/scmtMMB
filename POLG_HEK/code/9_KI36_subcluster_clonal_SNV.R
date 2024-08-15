setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code")
source("../../global_func/variant_calling.R")
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(Matrix)

mmat <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KI36.rds")
md <- read.csv("../output/1_HEK_metadata.csv") %>% 
  filter(condition == "KI36" & seurat_clusters %in% c(3,4,5,7)) # C2 = 5, C1 = 3+4, C3 =7
mitoSE <- readRDS("../POLG_HEK_large_data_files/input/mgatk/mgatk.rds")


# identify high-confident subcluster specific variants
barcode <- md %>% group_by(seurat_clusters) %>% summarize(barcode = list(X)) %>% tibble::deframe()
C1 <- call_mutations_mgatk(mitoSE[,c(barcode$`3`, barcode$`4`)])
C2 <- call_mutations_mgatk(mitoSE[,barcode$`5`])
C3 <- call_mutations_mgatk(mitoSE[,barcode$`7`])

mtMut <- function(mmat, n_cells = 1, vmr.thres = 0.01){
  var <- as.data.frame(rowData(mmat))
  var <- var %>% subset(n_cells_conf_detected >= n_cells &
                          strand_correlation > 0.65 &
                          vmr >= vmr.thres &
                          variant %ni% c("301A>C", "302A>C", "310T>C", "316G>C")) %>% pull(variant)
  mmat <- mmat[var,]
  return(mmat)
}

C1.f <- mtMut(C1, n_cells = 50)
C2.f <- mtMut(C2, n_cells = 50)
C3.f <- mtMut(C3, n_cells = 50)

Venn <- BioVenn::draw.venn(
  rownames(C1.f), rownames(C2.f), rownames(C3.f),
  title = NULL, subtitle = NULL,
  xtitle = NULL, ytitle = NULL, ztitle = NULL, 
  nr_c = "black", x_c = "#960096", y_c = "#960000", z_c = "#960055",
  bg_c = "white", output = "pdf", filename = "../plot/9_subcluster_KI36.pdf", width = 800, height =800
)

#MT-CO2
as.data.frame(rowData(C1.f)[Venn$x_only,]) %>% arrange(desc(mean)) %>% filter(position %in% 7586:8269)
as.data.frame(rowData(C2.f)[Venn$y_only,]) %>% arrange(desc(mean)) %>% filter(position %in% 7586:8269)
as.data.frame(rowData(C3.f)[Venn$z_only,]) %>% arrange(desc(mean)) %>% filter(position %in% 7586:8269)
#MT-CYB
as.data.frame(rowData(C1.f)[Venn$x_only,]) %>% arrange(desc(mean)) %>% filter(position %in% 14747:15887)
as.data.frame(rowData(C2.f)[Venn$y_only,]) %>% arrange(desc(mean)) %>% filter(position %in% 14747:15887)
as.data.frame(rowData(C3.f)[Venn$z_only,]) %>% arrange(desc(mean)) %>% filter(position %in% 14747:15887)
#MT-ND4L
as.data.frame(rowData(C1.f)[Venn$x_only,]) %>% arrange(desc(mean)) %>% filter(position %in% 10470:10766)
as.data.frame(rowData(C2.f)[Venn$y_only,]) %>% arrange(desc(mean)) %>% filter(position %in% 10470:10766)
as.data.frame(rowData(C3.f)[Venn$z_only,]) %>% arrange(desc(mean)) %>% filter(position %in% 10470:10766)



VAF.long <- as.data.frame(assay(mmat)) %>%
  # rename_with(.fn = ~paste0("Cell_", seq_along(.))) %>%
  rownames_to_column("variant") %>%
  separate(variant, into = c("Position", "Change"), sep = "(?<=\\d)(?=\\D)") %>%
  mutate(Position = as.numeric(Position)) %>%
  pivot_longer(cols = -c("Position", "Change"), names_to = "ID", values_to = "VAF")

VAF.long <- inner_join(VAF.long, md[,c("X","seurat_clusters", "condition", "atacUMAP_1", "atacUMAP_2")], by = c("ID" = "X"))

VAF.long <- VAF.long %>% mutate(present = VAF>0) %>% filter(atacUMAP_1 < 4 & atacUMAP_2 >1.9)

plot_mtSNV <- function(pos, var){
  
  p1 <- VAF.long %>% filter(Position == pos & Change == var & atacUMAP_1 < 4 & atacUMAP_2 >1.9) %>% arrange(VAF) %>% 
    ggplot(aes(atacUMAP_1, atacUMAP_2, color=VAF)) + geom_point() + scale_color_viridis_c(limits = c(0, 1), breaks = seq(0,1,0.2)) +
    ggtitle(paste0(pos, var)) + theme_void()
  
  p2 <- VAF.long %>% filter(Position == pos & Change == var & atacUMAP_1 < 4 & atacUMAP_2 >1.9) %>% arrange(present) %>% 
    ggplot(aes(atacUMAP_1, atacUMAP_2, color=present)) + geom_point() + scale_color_manual(values = c("grey", "red")) + 
    theme_void() + ggtitle(paste0(pos, var))
  
  p <- cowplot::plot_grid(p1, p2,ncol = 1)
  ggsave(plot = p, paste0("../plot/9_subcluster_clonal_var",pos, gsub(">", "", var),".pdf"), width = 5, height=5)
  p
  }

plot_mtSNV(7667, "C>T") # MT-CO2
plot_mtSNV(7927, "C>T") # MT-CO2
plot_mtSNV(15596, "G>A") # MT-CYB
plot_mtSNV(15378, "T>C") # MT-CYB
plot_mtSNV(10680, "G>A") # MT-ND4L
plot_mtSNV(10670, "C>T") # MT-ND4L
