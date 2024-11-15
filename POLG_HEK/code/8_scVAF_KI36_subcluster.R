setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code")
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(Signac)
library(Seurat)

KI36_clone <- read.csv("../output/scVAF/KI36.csv", row.names = "X")
mmtx <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KI36.rds")
md <- read.csv("../output/1_HEK_metadata.csv", row.names = "X") %>% 
  dplyr::filter(condition == "KI36" & seurat_clusters %in% c(3,4,5,7)) %>% # C2 = 5, C1 = 3+4, C3 =7
  mutate(subcluster = case_when(seurat_clusters == "5" ~ "C2",
                                seurat_clusters == "7" ~ "C3",
                                TRUE ~ "C1"),
         barcode = rownames(.))

md$T15378C <- assay(mmtx)["15378T>C",rownames(md)]>0.1
md$C7667T <- assay(mmtx)["7667C>T",rownames(md)]>0.1
md$C7927T <- assay(mmtx)["7927C>T",rownames(md)]>0.1

md$clonotype <- KI36_clone[rownames(md),"seurat_clusters"]



table(md$T15378C, md$subcluster)
#         C1   C2   C3
# FALSE 1071  413   36
# TRUE    12    7   78
table(md$C7667T, md$subcluster)
#         C1  C2  C3
# FALSE 185 403 105
# TRUE  898  17   9
table(md$C7927T, md$subcluster)
#         C1   C2   C3
# FALSE 1058   35  101
# TRUE    25  385   13

md_sub <- md %>% dplyr::filter(
  (T15378C == "TRUE" & md$subcluster == "C3") | 
    (C7927T == "TRUE" & md$subcluster == "C2") |
    (C7667T == "TRUE" & md$subcluster == "C1") ) %>% 
  group_by(subcluster) %>%
  slice_head(n=10) 


lapply(unique(md$subcluster), function(sc){
  
  CellID <- md_sub %>% 
    dplyr::filter(subcluster == sc) %>% 
    pull(barcode)
  
  VAF.long <- as.data.frame(assay(mmtx)[,CellID]) %>% 
    rename_with(.fn = ~paste0("Cell_", seq_along(.))) %>%
    rownames_to_column("variant") %>% 
    separate(variant, into = c("Position", "Change"), sep = "(?<=\\d)(?=\\D)") %>% 
    mutate(Position = as.numeric(Position)) %>% 
    pivot_longer(cols = -c("Position", "Change"), names_to = "ID", values_to = "VAF") 
  
  #data <- as.matrix(assay(mmtx)[,CellID])
  # Transpose data to have documents as rows and terms as columns
  #data_t <- t(data)
  
  # Compute cosine similarity matrix
  #cosine_similarity <- proxy::simil(data_t, method = "cosine")
  #print(cosine_similarity)
  
  lapply(unique(VAF.long$ID),function(cell){
    p <- VAF.long %>% dplyr::filter(ID == cell & VAF != 0) %>% 
      ggplot() +
      geom_segment(aes(x=Position, xend=Position, y=0, yend=VAF)) +
      geom_point(aes(x = Position, y = VAF), size=0.1, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=2) +
      scale_x_continuous(limits = c(0,16569), expand = c(0,0)) +
      scale_y_continuous(limits = c(-0.5,1), expand = c(0,0)) +
      geom_hline(yintercept = c(0,1)) +
      geom_hline(yintercept = seq(0.1,0.9,0.1), linetype = "dashed",size= 0.2, color = "grey30") +
      geom_rect(data = mitochondrial_data, 
                aes(xmin = Starting, xmax = Ending, ymin = -0.2, ymax = -0.1, fill = color), alpha = 1) +
      scale_fill_manual(values = mtGene_colors) +
      coord_polar() + 
      guides(fill=guide_legend(ncol=2))+
      theme_void() + NoLegend() + theme(aspect.ratio=1/1) + 
      ggtitle(paste0(sc, "_",cell))
    ggsave(plot = p, filename = paste0("../plot/scVAF/KI36_subcluster/", sc, "_", cell,".pdf"), width = 5, height = 5)
  })
  
})



table(md$T15378C, md$clonotype) #C3
#         0   1   2   3   4   5   6   7   8
# FALSE 380 377 346 220   2  38  22  22  21
# TRUE   33   6  37  19  88   0   5   1   0

table(md$C7667T, md$clonotype) #mainly C1
#         0   1   2   3   4   5   6   7   8
# FALSE  20 361  25  18  88   7   3  23  19
# TRUE  393  22 358 221   2  31  24   0   2

table(md$C7927T, md$clonotype) #mainly C2
#         0   1   2   3   4   5   6   7   8
# FALSE 397   5 371 227  88  34   0   7  20
# TRUE   16 378  12  12   2   4  27  16   1

