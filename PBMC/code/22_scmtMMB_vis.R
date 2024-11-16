setwd("~/Ludwig_lab/scmtMMB/PBMC/code/")
library(ggplot2)
library(scales)
library(dplyr)
'%ni%' <- Negate('%in%')


scmtMMB <- read.csv("../output/2_scmtMMB.csv", row.names = "X")
scmtMMB$X3243 <- "+"
scmtMMB_ex <- read.csv("../output/2_scmtMMB_exclude_3243.csv", row.names = "X") 
scmtMMB_ex$X3243 <- "-"

color_vec <- c("M80" = "darkorchid4","M60" = "firebrick4", "M35" = "violetred3", 
               "M29" = "firebrick1", "H47" = "blue3", "H05" = "dodgerblue3")

P1 <- rbind(scmtMMB %>% filter(!sample %in% c("H05", "H47")), scmtMMB_ex) %>% dplyr::filter(symbol == "Total") %>%
  mutate(X3243 = factor(X3243, level = c("+", "-"))) %>% 
  group_by(sample, X3243) %>%
  mutate(
    q5 = quantile(mutation_per_MB, 0.05),
    q95 = quantile(mutation_per_MB, 0.95)
  ) %>%
  filter(mutation_per_MB >= q5, mutation_per_MB <= q95) %>%  # Filter out values outside 5th and 95th percentiles
  ungroup() %>% 
  ggplot(aes(x = sample, y = mutation_per_MB)) +
  #geom_jitter(aes(fill = sample), size = 0.1, color = "grey") + 
  geom_violin(aes(fill = X3243), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = c("firebrick", "darkslategray3")) +
  # geom_boxplot(aes(fill = X3243), scale = "width", alpha = 0.8, width = 0.1, 
  #             position = position_dodge(width=0.9), outlier.size = 0.1) +
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), 
             aes(group = interaction(sample, X3243)), shape = 23, size = 1, fill = "white",
             position = position_dodge(width=0.9)) + # overlay bulk data as point
  theme_bw() + xlab(NULL) + ylab("Total scMPM")+
  theme(aspect.ratio=1/1, axis.title.x=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
P1.nl <- P1 + Seurat::NoLegend() 
P1.leg <- ggpubr::get_legend(P1) %>% ggpubr::as_ggplot()

ggsave(plot = P1.nl, "../plot/22_scMPM_total_pm_MELAS_q5_q95.pdf", width=3, height=3)
ggsave(plot = P1.leg, "../plot/22_scMPM_total_pm_MELAS_leg_q5_q95.pdf", width=3, height=3)


P2 <- rbind(scmtMMB %>% filter(!sample %in% c("H05", "H47")), scmtMMB_ex) %>% dplyr::filter(symbol == "Total") %>%
  mutate(X3243 = factor(X3243, level = c("+", "-"))) %>% 
  group_by(sample, X3243) %>%
  mutate(
    q5 = quantile(MSS_weighted, 0.05),
    q95 = quantile(MSS_weighted, 0.95)
  ) %>%
  filter(MSS_weighted >= q5, MSS_weighted <= q95) %>%  # Filter out values outside 5th and 95th percentiles
  ungroup() %>% 
  ggplot(aes(x = sample, y = MSS_weighted)) +
  #geom_jitter(aes(fill = sample), size = 0.1, color = "grey") + 
  geom_violin(aes(fill = X3243), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = c("firebrick", "darkslategray3")) +
  # geom_boxplot(aes(fill = X3243), scale = "width", alpha = 0.8, width = 0.1, 
  #             position = position_dodge(width=0.9), outlier.size = 0.1) +
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), 
             aes(group = interaction(sample, X3243)), shape = 23, size = 1, fill = "white",
             position = position_dodge(width=0.9)) + # overlay bulk data as point
  theme_bw() + xlab(NULL) + ylab("Total scwMSS")+
  theme(aspect.ratio=1/1, axis.title.x=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
P2.nl <- P2 + Seurat::NoLegend() 
P2.leg <- ggpubr::get_legend(P2) %>% ggpubr::as_ggplot()

ggsave(plot = P2.nl, "../plot/22_scwMSS_total_pm_MELAS_q5_q95.pdf", width=3, height=3)



# Visulization
theme_1_nogrid <- theme(aspect.ratio=1/1, axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

### Weighted MSS
p11 <-scmtMMB %>% dplyr::filter(symbol == "Total") %>%
  ggplot(aes(x = sample, y = MSS_weighted)) +
  #geom_jitter(aes(fill = sample), size = 0.1, color = "grey") + 
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  theme_bw() + xlab(NULL) + ylab("Total scwMSS")+
  theme_1_nogrid + Seurat::NoLegend()
# facet_wrap(~factor(symbol, levels = c("CI", "CIII", "CIV", "CV", "tRNA", "rRNA", "Total")), ncol =3) + 
#scale_y_continuous(trans = log2_trans(), 
#                   breaks = trans_breaks("log2", function(x) 2^x), 
#                   labels = trans_format("log2", math_format(2^.x))) + 
# labs(y = bquote(paste(log[2], "(Weighted MSS + 0.01)")))
p11
ggsave(plot = p11, "../plot/22_scwMSS_total_MELAS.pdf", width=3, height=3)


p112 <- scmtMMB %>% dplyr::filter(symbol %in% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA")) %>%
  ggplot(aes(x = sample, y = MSS_weighted + 0.01)) +
  #geom_jitter(aes(fill = sample), size = 0.1, color = "grey") + 
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  facet_wrap(~ symbol) + # overlay bulk data as point
  theme_bw() + xlab(NULL) + theme_1_nogrid +
scale_y_continuous(trans = log2_trans(), 
                  breaks = trans_breaks("log2", function(x) 2^x), 
                  labels = trans_format("log2", math_format(2^.x))) + 
labs(y = bquote(paste(log[2], "(Weighted MSS + 0.01)"))) + Seurat::NoLegend()
p112
ggsave(plot = p112, "../plot/22_scwMSS_complex_MELAS.pdf", width=6, height=4)

#exclude 3243A>G
p12 <- scmtMMB_ex %>% dplyr::filter(symbol == "Total") %>%
  ggplot(aes(x = sample, y = MSS_weighted)) +
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  theme_bw() + xlab(NULL) + ylab("Total scwMSS (exclude 3243A>G)")+
  theme_1_nogrid + Seurat::NoLegend()
p12
ggsave(plot = p12, "../plot/22_scwMSS_total_MELAS_ex3243.pdf", width=3, height=3)

p122 <- scmtMMB_ex %>% dplyr::filter(symbol == "tRNA") %>% # %>% dplyr::filter(symbol %in% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA")) %>%
  ggplot(aes(x = sample, y = MSS_weighted + 0.01)) +
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  facet_wrap(~ symbol) + # overlay bulk data as point
  theme_bw() + xlab(NULL) + theme_1_nogrid +
  scale_y_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x), 
                     labels = trans_format("log2", math_format(2^.x))) + 
  labs(y = bquote(paste(log[2], "(Weighted MSS + 0.01)"))) + Seurat::NoLegend()
p122
ggsave(plot = p122, "../plot/22_scwMSS_complex_MELAS_ex3243.pdf", width=3, height=3)

p13 <- scmtMMB_ex %>% dplyr::filter(symbol %ni% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA", "Total") & !grepl("MT_T",symbol)) %>%
  ggplot(aes(x = sample, y = MSS_weighted + 0.01)) +
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  facet_wrap(~ symbol, ncol =5) + # overlay bulk data as point
  theme_bw() + xlab(NULL) + theme_1_nogrid +
  scale_y_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x), 
                     labels = trans_format("log2", math_format(2^.x))) + 
  labs(y = bquote(paste(log[2], "(Weighted MSS + 0.01)"))) + Seurat::NoLegend()
p13
ggsave(plot = p13, "../plot/23_scwMSS_gene_MELAS_ex3243.pdf", width=10, height=6)

scmtMMB_ex %>% dplyr::filter(symbol %in% c("MT_ND3", "MT_ND4L")) %>%
  ggplot(aes(x = sample, y = MSS_weighted + 0.01)) +
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  facet_wrap(~ symbol, ncol =5) + # overlay bulk data as point
  theme_bw() + xlab(NULL) + theme_1_nogrid +
  scale_y_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x), 
                     labels = trans_format("log2", math_format(2^.x))) + 
  labs(y = bquote(paste(log[2], "(Weighted MSS + 0.01)"))) + Seurat::NoLegend()


## Mutations per MB
p21 <- scmtMMB %>% dplyr::filter(symbol == "Total") %>%
  ggplot(aes(x = sample, y = mutation_per_MB)) +
  #geom_jitter(aes(fill = sample), size = 0.1, color = "grey") + 
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  theme_bw() + xlab(NULL) + ylab("Total scMPM")+
  theme_1_nogrid+ Seurat::NoLegend() 
#scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
#                   labels = trans_format("log10", math_format(10^.x))) + 
#labs(y = bquote(paste(log[10], "(Mutation per MB + 1)")))

ggsave(plot = p21, "../plot/22_scMMB_total_MELAS.pdf", width=3, height=3)


P212 <- scmtMMB %>% dplyr::filter(symbol %in% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA")) %>%
  ggplot(aes(x = sample, y = mutation_per_MB + 1)) +
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  facet_wrap(~ symbol) +
  theme_bw() + xlab(NULL) + #ylab("Total scMMB (exclude 3243A>G)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) + 
  labs(y = bquote(paste(log[10], "(Mutation per MB + 1)"))) +
  theme_1_nogrid + Seurat::NoLegend()
P212
ggsave(plot = P212, "../plot/22_scMMB_complex_MELAS.pdf", width=6, height=4)


p22 <- scmtMMB_ex %>% dplyr::filter(symbol == "Total") %>%
  ggplot(aes(x = sample, y = mutation_per_MB)) +
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  theme_bw() + xlab(NULL) + ylab("Total scMMB (exclude 3243A>G)") +
  theme_1_nogrid+ Seurat::NoLegend() 
ggsave(plot = p22, "../plot/22_scMMB_total_MELAS_ex3243.pdf", width=3, height=3)



p222 <- scmtMMB_ex %>% dplyr::filter(symbol == "tRNA") %>% #%>% dplyr::filter(symbol %in% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA")) %>%
  ggplot(aes(x = sample, y = mutation_per_MB + 1)) +
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  facet_wrap(~ symbol) +
  theme_bw() + xlab(NULL) + #ylab("Total scMMB (exclude 3243A>G)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) + 
  labs(y = bquote(paste(log[10], "(Mutation per MB + 1)"))) +
  theme_1_nogrid+ Seurat::NoLegend() 

ggsave(plot = p222, "../plot/22_scMMB_complex_MELAS_ex3243.pdf", width=3, height=3)

p23 <- scmtMMB_ex %>% dplyr::filter(symbol %ni% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA", "Total") & !grepl("MT_T",symbol)) %>%
  ggplot(aes(x = sample, y = mutation_per_MB + 1)) +
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  facet_wrap(~ symbol, ncol = 5) +
  theme_bw() + xlab(NULL) + #ylab("Total scMMB (exclude 3243A>G)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) + 
  labs(y = bquote(paste(log[10], "(Mutation per MB + 1)"))) +
  theme_1_nogrid+ Seurat::NoLegend() 

ggsave(plot = p23, "../plot/23_scMMB_gene_MELAS_ex3243.pdf", width=10, height=6)

scmtMMB_ex %>% dplyr::filter(symbol %in% c("MT_ND4L", "MT_ND3")) %>%
  ggplot(aes(x = sample, y = mutation_per_MB + 1)) +
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  facet_wrap(~ symbol, ncol = 5) +
  theme_bw() + xlab(NULL) + #ylab("Total scMMB (exclude 3243A>G)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) + 
  labs(y = bquote(paste(log[10], "(Mutation per MB + 1)"))) +
  theme_1_nogrid+ Seurat::NoLegend() 


scmtMMB %>% dplyr::filter(symbol %in% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA")) %>%
  ggplot(aes(x = mutation_per_MB, color = sample)) +
  stat_ecdf() + scale_color_manual(values = color_vec) + facet_wrap(~ symbol, scales = "free_x") +
  theme_classic()  + Seurat::NoLegend()

# project scmtMMB on UMAP space:at functional group level

#### exclude 3243A>G
lapply(c("CI", "CIII", "CIV", "CV", "tRNA", "rRNA", "Total"), function(cat){
  
  df <- scmtMMB_ex %>% dplyr::filter(symbol == cat & barcode != "pseudobulk") #%>% 
  # mutate(log10_MpMB = log10(mutation_per_MB), log2_MSSw = log2(MSS_weighted)) %>% 
  # mutate(log10_MpMB = ifelse(is.infinite(log10_MpMB), NA, log10_MpMB),
  #        log2_MSSw = ifelse(is.infinite(log2_MSSw), NA, log2_MSSw)) %>% arrange(mutation_per_MB)
  
  p1 <- ggplot(df %>% arrange(mutation_per_MB), aes(x = refUMAP_1, y = refUMAP_2)) +
    geom_point(color = "grey", size = 0.01) +
    geom_point(aes(color = mutation_per_MB),size = 0.01) +
    scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra"),
                          limit = quantile(df$mutation_per_MB, c(0.02,0.98), na.rm =TRUE), na.value = "transparent") + 
    facet_wrap(~sample, nrow = 3) +
    theme(aspect.ratio=1) # +  ggtitle(paste0("MpMB_", cat)) 
  p1.nl <- p1 + theme_void() + Seurat::NoLegend() 
  p1.leg <- ggpubr::get_legend(p1) %>% ggpubr::as_ggplot()
  ggsave(plot = p1.nl, paste0("../plot/2_MpMB_UMAP_ex3243/MpMB_", cat,".pdf"), width = 4, height = 6)
  ggsave(plot = p1.leg, paste0("../plot/2_MpMB_UMAP_ex3243/MpMB_", cat,".leg.pdf"), width = 3, height = 2)
  
  
  p2 <- ggplot(df %>% arrange(MSS_weighted), aes(x = refUMAP_1, y = refUMAP_2)) +
    geom_point(color = "grey", size = 0.01) + 
    geom_point(aes(color = MSS_weighted), size = 0.01) + 
    scale_color_viridis_c(limit = quantile(df$MSS_weighted, c(0.02,0.98), na.rm =TRUE), 
                          na.value = "transparent") +
    theme(aspect.ratio=1) + facet_wrap(~sample, nrow = 3) 
  
  p2.nl <- p2 + theme_void() + Seurat::NoLegend()
  p2.leg <- ggpubr::get_legend(p2) %>% ggpubr::as_ggplot()
  ggsave(plot = p2.nl, paste0("../plot/2_MSSw_UMAP_ex3243/MSSw_", cat,".pdf"), width = 4, height = 6)
  ggsave(plot = p2.leg, paste0("../plot/2_MSSw_UMAP_ex3243/MSSw_", cat,".leg.pdf"), width = 3, height = 2)
  
}) 

#### include 3243A>G

lapply(c("CI", "CIII", "CIV", "CV", "tRNA", "rRNA", "Total"), function(cat){
  
  df <- scmtMMB %>% dplyr::filter(symbol == cat & barcode != "pseudobulk") #%>% 
  # mutate(log10_MpMB = log10(mutation_per_MB), log2_MSSw = log2(MSS_weighted)) %>% 
  # mutate(log10_MpMB = ifelse(is.infinite(log10_MpMB), NA, log10_MpMB),
  #        log2_MSSw = ifelse(is.infinite(log2_MSSw), NA, log2_MSSw)) %>% arrange(mutation_per_MB)
  
  p1 <- ggplot(df %>% arrange(mutation_per_MB), aes(x = refUMAP_1, y = refUMAP_2)) +
    geom_point(color = "grey", size = 0.01) +
    geom_point(aes(color = mutation_per_MB),size = 0.01) +
    scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra"),
                          limit = quantile(df$mutation_per_MB, c(0.02,0.98), na.rm =TRUE), na.value = "transparent") + 
    facet_wrap(~sample, nrow = 3) +
    theme(aspect.ratio=1) # +  ggtitle(paste0("MpMB_", cat)) 
  p1.nl <- p1 + theme_void() + Seurat::NoLegend() 
  p1.leg <- ggpubr::get_legend(p1) %>% ggpubr::as_ggplot()
  ggsave(plot = p1.nl, paste0("../plot/2_MpMB_UMAP/MpMB_", cat,".pdf"), width = 4, height = 6)
  ggsave(plot = p1.leg, paste0("../plot/2_MpMB_UMAP/MpMB_", cat,".leg.pdf"), width = 3, height = 2)
  
  
  p2 <- ggplot(df %>% arrange(MSS_weighted), aes(x = refUMAP_1, y = refUMAP_2)) +
    geom_point(color = "grey", size = 0.01) + 
    geom_point(aes(color = MSS_weighted), size = 0.01) + 
    scale_color_viridis_c(limit = quantile(df$MSS_weighted, c(0.02,0.98), na.rm =TRUE), 
                          na.value = "transparent") +
    theme(aspect.ratio=1) + facet_wrap(~sample, nrow = 3) 
  
  p2.nl <- p2 + theme_void() + Seurat::NoLegend()
  p2.leg <- ggpubr::get_legend(p2) %>% ggpubr::as_ggplot()
  ggsave(plot = p2.nl, paste0("../plot/2_MSSw_UMAP/MSSw_", cat,".pdf"), width = 4, height = 6)
  ggsave(plot = p2.leg, paste0("../plot/2_MSSw_UMAP/MSSw_", cat,".leg.pdf"), width = 3, height = 2)
  
}) 

M29_3243 <- scmtMMB %>% filter(sample == "M29", symbol == "Total") %>% 
  arrange(X3243A_G) %>% 
  ggplot(aes(x=refUMAP_1, y = refUMAP_2, col = X3243A_G)) +
  geom_point(size = 0.01) + viridis::scale_color_viridis(begin = 0, end = 1)+
  theme_classic() + theme_void() + Seurat::NoLegend()
ggsave(plot = M29_3243, paste0("../plot/2_M29_3243.pdf"), width = 2, height = 2)

M35_3243 <- scmtMMB %>% filter(sample == "M35", symbol == "Total") %>% 
  arrange(X3243A_G) %>% 
  ggplot(aes(x=refUMAP_1, y = refUMAP_2, col = X3243A_G)) +
  geom_point(size = 0.01) + viridis::scale_color_viridis(begin = 0, end = 1)+
  theme_classic() + theme_void() + Seurat::NoLegend()
ggsave(plot = M35_3243, paste0("../plot/2_M35_3243.pdf"), width = 2, height = 2)


# scmtMMB %>% filter(grepl("M", sample) & symbol == "Total") %>% 
#   ggplot(aes(x = X3243A_G, color = sample)) +
#   stat_ecdf() + scale_color_manual(values = color_vec) + 
#   scale_x_continuous(limits = c(-0.05, 1.01), expand = c(0,0)) +
#   scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
#   theme_classic() + theme_cor + Seurat::NoLegend() +
#   xlab("m.3243A>G Heteroplasmy") +
#   ylab("Cumulative frequency")

# scmtMMB_ex %>% dplyr::filter(symbol %in% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA")) %>%
#   ggplot(aes(x = sample, y = MSS_weighted + 0.01)) +
#   geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
#   scale_fill_manual(values = color_vec) +
#   geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), 
#              aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
#   facet_wrap(~ symbol, scales = "free_y") +
#   theme_bw() + xlab(NULL) + ylab("Total scwMSS (exclude 3243A>G)")+
#   theme_1_nogrid
