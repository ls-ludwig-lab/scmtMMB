setwd("~/Ludwig_lab/scmtMMB/PBMC/code/")
library(Matrix)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(viridis)
library(BuenColors)
library(ggrepel)

scmtMMB <- read.csv("../output/2_scmtMMB.csv", row.names = "X") %>% filter(barcode != "pseudobulk")

scmtMMB <- read.csv("../output/2_scmtMMB.csv", row.names = "X") %>% 
  filter(symbol == "Total" & barcode != "pseudobulk") %>%
  mutate(
    celltype = case_when(
      predicted.celltype.l2 %in% c("CD8 TEM", "CD8 TCM") ~ "CD8 Eff/Mem",
      predicted.celltype.l2 %in% c("CD4 TEM", "CD4 TCM", "CD4 CTL") ~ "CD4 Eff/Mem",
      predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory", "Plasmablast") ~ "B cell",
      predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC") ~ "Monocyte/DC",
      grepl("NK|ILC", predicted.celltype.l2) ~ "NK/ILC",
      predicted.celltype.l2 %in% c("Doublet", "Eryth", "Platelet", "ASDC", NA) ~ "discard",
      TRUE ~ predicted.celltype.l2
    )
  ) %>% filter(celltype != "discard")


##########  H05  ###############

H5_mmat <- readRDS("../PBMC_large_data_files/output/01_mmat_H05.rds")

pos_ND4L <- 10470:10766
idx <- rowData(H5_mmat)$position %in% pos_ND4L
H5_mmat_ND4L <- H5_mmat[idx,]
p1 <- as.data.frame(rowData(H5_mmat_ND4L)) %>% arrange(desc(mean)) %>% 
  mutate(rank = 1:25) %>% ggplot(aes(x = rank, y = mean)) + geom_point(size = 0.5) + theme_classic() +
  geom_text_repel(aes(label = variant)) +
  ggtitle("Mean VAF of MT-ND4L variants in H05") #10599G>A
p1
ggsave(plot = p1, paste0("../plot/25_H05_ND4L.pdf"), width = 2, height = 2)

df_H05 <- scmtMMB %>% filter(sample == "H05" & symbol == "Total") %>% 
  mutate(VAF = assay(H5_mmat)["10599G>A",.$barcode]) %>% arrange(VAF)

p2 <- df_H05 %>% 
  ggplot(aes(x=refUMAP_1, y = refUMAP_2, col = VAF)) +
  geom_point(size = 0.01) + scale_color_viridis(begin = 0, end = 1)+
  theme_classic() + theme_void() + Seurat::NoLegend()
p2

ggsave(plot = p2, paste0("../plot/25_H05_ND4L_G10559A.pdf"), width = 2, height = 2)

# density plot
df_H05$af_weight <- Nebulosa:::calculate_density(df_H05$VAF, cbind(df_H05$refUMAP_1, df_H05$refUMAP_2),
                                                  method = "wkde", adjust = 5)
p2.2 <- ggplot(df_H05 %>% arrange(af_weight),
                               aes(x = refUMAP_1, y = refUMAP_2, color = af_weight)) +
  geom_point(size = 0.1) +
  scale_color_gradientn(colors = BuenColors::jdb_palette("flame_flame")) + theme_void() 

p2.2.nl <- p2.2 + theme(legend.position = "none")
p2.2.leg <- ggpubr::get_legend(p2.2) %>% ggpubr::as_ggplot()

ggsave(plot = p2.2.leg, paste0("../plot/25_density.leg.pdf"), width = 3, height = 2)

ggsave(plot = p2.2.nl, paste0("../plot/25_H05_ND4L_G10559A_smooth_umap.pdf"), width = 2, height = 2)

##########  M60  ###############

M60_mmat <- readRDS("../PBMC_large_data_files/output/01_mmat_M60.rds")

pos_ND3 <- 10059:10404
idx <- rowData(M60_mmat)$position %in% pos_ND3
M60_mmat_ND3 <- M60_mmat[idx,]

p3 <-as.data.frame(rowData(M60_mmat_ND3)) %>% arrange(desc(mean)) %>% 
  mutate(rank = 1:nrow(.)) %>% ggplot(aes(x = rank, y = mean)) + geom_point(size =0.5) + theme_classic() +
  geom_text_repel(aes(label = variant)) +
  ggtitle("Mean VAF of MT-ND3 variants in M60") #10270T>C 
p3
ggsave(plot = p3, paste0("../plot/25_M60_ND3.pdf"), width = 2, height = 2)

df_M60 <- scmtMMB %>% filter(sample == "M60" & symbol == "Total") %>% 
  mutate(VAF = assay(M60_mmat_ND3)["10270T>C",.$barcode]) %>% arrange(VAF)
  
p4 <- df_M60 %>% 
  ggplot(aes(x=refUMAP_1, y = refUMAP_2, col = VAF)) +
  geom_point(size = 0.01) + scale_color_viridis(begin = 0, end = 1)+
  theme_classic()+ theme_void() + Seurat::NoLegend()
p4
ggsave(plot = p4, paste0("../plot/25_M60_ND3_T10270C.pdf"), width = 2, height = 2)

# density plot
df_M60$af_weight <- Nebulosa:::calculate_density(df_M60$VAF, cbind(df_M60$refUMAP_1, df_M60$refUMAP_2),
                                                 method = "wkde", adjust = 5)
p2.2 <- ggplot(df_M60 %>% arrange(af_weight),
               aes(x = refUMAP_1, y = refUMAP_2, color = af_weight)) +
  geom_point(size = 0.1) +
  scale_color_gradientn(colors = BuenColors::jdb_palette("flame_flame")) + theme_void() +
  theme(legend.position = "none")
p2.2
ggsave(plot = p2.2, paste0("../plot/25_M60_ND3_T10270C_smooth_umap.pdf"), width = 2, height = 2)


##########  M80  ###############
M80_mmat <- readRDS("../PBMC_large_data_files/output/01_mmat_M80.rds")

idx <- rowData(M80_mmat)$position %in% pos_ND3
M80_mmat_ND3 <- M80_mmat[idx,]
p5 <-  as.data.frame(rowData(M80_mmat_ND3)) %>% arrange(desc(mean)) %>% 
  mutate(rank = 1:nrow(.)) %>% ggplot(aes(x = rank, y = mean)) + 
  geom_point(size =0.5) + geom_text_repel(aes(label = variant)) +
  theme_classic() +
  ggtitle("Mean VAF of MT-ND3 variants in M80") # 10398A>G
p5 
ggsave(plot = p5 , paste0("../plot/25_M80_ND3.pdf"), width = 2, height = 2)

df_M80 <- scmtMMB %>% filter(sample == "M80" & symbol == "Total") %>% 
  mutate(VAF = assay(M80_mmat_ND3)["10398A>G",.$barcode]) %>% arrange(VAF)
  
p6 <- df_M80 %>% 
  ggplot(aes(x=refUMAP_1, y = refUMAP_2, col = VAF)) +
  geom_point(size = 0.01) + scale_color_viridis(begin = 0, end = 1)+
  theme_classic() + theme_void() + Seurat::NoLegend()
p6
ggsave(plot = p6, paste0("../plot/25_M80_ND3_A10398G.pdf"), width = 2, height = 2)

# density plot
df_M80$af_weight <- Nebulosa:::calculate_density(df_M80$VAF, cbind(df_M80$refUMAP_1, df_M80$refUMAP_2),
                                                 method = "wkde", adjust = 5)
p2.2 <- ggplot(df_M80 %>% arrange(af_weight),
               aes(x = refUMAP_1, y = refUMAP_2, color = af_weight)) +
  geom_point(size = 0.1) +
  scale_color_gradientn(colors = BuenColors::jdb_palette("flame_flame")) + theme_void() +
  theme(legend.position = "none")
p2.2
ggsave(plot = p2.2, paste0("../plot/25_M80_ND3_A10398G_smooth_umap.pdf"), width = 2, height = 2)


scmtMMB %>% filter(sample == "M80" & symbol == "MT_ND3") %>% 
  mutate(A10398G = assay(M80_mmat)["10398A>G",.$barcode]) %>% 
  arrange(desc(mtDNA_depth)) %>% 
  ggplot(aes(x=A10398G, y = mutation_per_MB)) +
  geom_point(size = 0.01) +# scale_color_viridis() +
  theme_classic() 

scmtMMB_ex %>% filter(sample == "M80" & symbol == "MT_ND3") %>% 
  mutate(A10398G = assay(M80_mmat)["10398A>G",.$barcode]) %>% 
  arrange(desc(mtDNA_depth)) %>% 
  ggplot(aes(x=A10398G, y = MSS_weighted)) +
  geom_point(size = 0.01) +# scale_color_viridis() +
  theme_classic() 


######################

pos_RNR2 <- 1671:3229
idx <- rowData(M80_mmat)$position %in% pos_RNR2
M80_mmat_RNR2 <- M80_mmat[idx,]
p7 <- as.data.frame(rowData(M80_mmat_RNR2)) %>% arrange(desc(mean)) %>% 
  mutate(rank = 1:nrow(.)) %>% ggplot(aes(x = rank, y = mean)) + geom_point(size = 0.5) + theme_classic() +
  geom_text_repel(aes(label = variant)) +
  ggtitle("Mean VAF of MT-RNR2 variants in M80") #1969G>A, 3079G>A #2702G>A
p7
ggsave(plot = p7, paste0("../plot/25_M80_RNR2.pdf"), width = 2, height = 2)

p8.1 <- scmtMMB %>% filter(sample == "M80" & symbol == "Total") %>% 
  mutate(VAF = assay(M80_mmat_RNR2)["1969G>A",.$barcode]) %>% arrange(VAF) %>% 
  ggplot(aes(x=refUMAP_1, y = refUMAP_2, col = VAF)) +
  geom_point(size = 0.01) + scale_color_viridis(begin = 0, end = 1)+
  theme_classic()+ theme_void() + Seurat::NoLegend()
p8.1
ggsave(plot = p8.1, paste0("../plot/25_M80_RNR2_G1969A.pdf"), width = 2, height = 2)

p8.2 <- scmtMMB %>% filter(sample == "M80" & symbol == "Total") %>% 
  mutate(VAF = assay(M80_mmat_RNR2)["3079G>A",.$barcode]) %>% arrange(VAF) %>% 
  ggplot(aes(x=refUMAP_1, y = refUMAP_2, col = VAF)) +
  geom_point(size = 0.01) + scale_color_viridis(begin = 0, end = 1)+
  theme_classic()+ theme_void() + Seurat::NoLegend()
p8.2
ggsave(plot = p8.2, paste0("../plot/25_M80_RNR2_G3079A.pdf"), width = 2, height = 2)

p8.3 <- scmtMMB %>% filter(sample == "M80" & symbol == "Total") %>% 
  mutate(VAF = assay(M80_mmat_RNR2)["2702G>A",.$barcode]) %>% arrange(VAF) %>% 
  ggplot(aes(x=refUMAP_1, y = refUMAP_2, col = VAF)) +
  geom_point(size = 0.01) + scale_color_viridis(begin = 0, end = 1)+
  theme_classic()+ theme_void() + Seurat::NoLegend()
p8.3
ggsave(plot = p8.3, paste0("../plot/3_M80_RNR2_G2702A.pdf"), width = 2, height = 2)