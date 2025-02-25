setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
source("../../global_func/quantifyMMB.R")
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(data.table)
'%ni%' <- Negate('%in%')


scmtMMB_HEK <- read.csv("../output/5_scmtMMB_HEK.csv", row.names = "X")
color_vec <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")

scmtMMB_HEK_com <- scmtMMB_HEK %>% filter(symbol %in% c("CI", "CIII", "CIV", "CV", "tRNA", "rRNA", "Total")) %>% 
  mutate(symbol = factor(symbol,levels = c("CI", "CIII", "CIV", "CV", "tRNA", "rRNA", "Total")))

p2 <- scmtMMB_HEK_com %>% 
  ggplot(aes(x = sample, y =  MSS_weighted + 0.01)) +
  geom_jitter(aes(fill = sample), size = 0.1, color = "grey") + 
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = scmtMMB_HEK_com %>% filter(barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 2, fill = "white") + # overlay bulk data as point
  facet_wrap(~symbol, ncol =10) + 
  theme_bw() + xlab(NULL) +
  scale_y_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x), 
                     labels = trans_format("log2", math_format(2^.x))) + 
  theme(aspect.ratio=1/1, axis.title.x=element_blank(), 
        axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(y = bquote(paste(log[2], "(scwMSS + 0.01)")))


ggsave(plot = p2, "../plot/52_MSS_weighted_complex.pdf", width = 14, height = 2.5)


p3 <- scmtMMB_HEK_com %>%
  ggplot(aes(x = sample, y = mutation_per_MB + 1)) +
  geom_jitter(aes(fill = sample), size = 0.1, color = "grey") + 
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% filter(barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 2, fill = "white") + # overlay bulk data as point
  facet_wrap(~symbol, ncol =10) + 
  theme_bw() + xlab(NULL) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) + 
  theme(aspect.ratio=1/1, axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(y = bquote(paste(log[10], "(scmtMPM + 1)")))

# check if the elevated CI/CIV MMB cells in CTRL are the same ones
scmtMMB_HEK_com %>% filter(sample == "CTRL" & symbol %in% c("CI","CIV")) %>% 
  select(barcode, symbol, mutation_per_MB) %>% 
  pivot_wider(names_from = symbol, values_from = mutation_per_MB) %>% 
  ggplot(aes(CI, CIV)) +
  geom_point() + theme_bw()

p4 <- scmtMMB_HEK %>% filter(symbol == "Total") %>% 
  ggplot(aes(x = log10(mtDNA_depth), y = mutation_per_MB + 1)) +
  geom_point(aes(col = sample)) + facet_wrap(~sample) + 
  scale_color_manual(values = color_vec) + 
  theme_bw() + ylab("Mutation per MB (Total)") + theme(legend.position="none") +
  geom_smooth(method = "lm"))
ggsave(plot = p4, paste0("../plot/54_CorDepth_MpMB.pdf")

p5 <- scmtMMB_HEK %>% filter(symbol == "Total") %>% ggplot(aes(x = log10(mtDNA_depth), y = MSS_weighted, col = sample)) +
      geom_point() + facet_wrap(~sample) + scale_color_manual(values = color_vec) + 
      theme_bw() + ylab("Weighted MSS (Total)") + theme(legend.position="none") +
      geom_smooth(method = "lm")
ggsave(plot = p5, paste0("../plot/55_CorDepth_MSSw.pdf"))

(p6 <- scmtMMB_HEK %>% filter(symbol == "Total") %>% 
    ggplot(aes(x= mutation_per_MB +1, y = MSS_weighted, color = sample)) +
    geom_point(size = 0.01) + scale_color_manual(values = color_vec) +
    facet_wrap(~sample, scales = "free", ncol = 3) + 
    theme_bw() + theme(legend.position="none") +
    geom_smooth(method = "lm")) 
ggsave(plot = p6, paste0("../plot/56_Cor_MpMB_MSSw.pdf"))

library(purrr)
scmtMMB_HEK %>% filter(symbol == "Total") %>% 
  group_by(sample) %>% 
  nest() %>% 
  mutate(
    #mod1 = map(data, ~ lm(mutation_per_MB ~ log10(mtDNA_depth), data = .x)),
    rp.MPM = map_dbl(data, ~ cor(.x$mutation_per_MB, log10(.x$mtDNA_depth), method = "pearson", use = "complete.obs")),
    rs.MPM = map_dbl(data, ~ cor(.x$mutation_per_MB, log10(.x$mtDNA_depth), method = "spearman", use = "complete.obs")),
    #mod2 = map(data, ~ lm(MSS_weighted ~ log10(mtDNA_depth), data = .x)),
    rp.MSS = map_dbl(data, ~ cor(.x$MSS_weighted, log10(.x$mtDNA_depth), method = "pearson", use = "complete.obs")),
    rs.MSS = map_dbl(data, ~ cor(.x$MSS_weighted, log10(.x$mtDNA_depth), method = "spearman", use = "complete.obs")),
    #mod3 = map(data, ~ lm(MSS_weighted ~ mutation_per_MB, data = .x)),
    rp.MMB = map_dbl(data, ~ cor(.x$MSS_weighted, .x$mutation_per_MB, method = "pearson", use = "complete.obs")),
    rs.MMB = map_dbl(data, ~ cor(.x$MSS_weighted, .x$mutation_per_MB, method = "spearman", use = "complete.obs"))
  ) %>% 
  ungroup()
# # A tibble: 3 × 8
#   sample data                  rp.MPM rs.MPM rp.MSS rs.MSS rp.MMB rs.MMB
#   <chr>  <list>                 <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
# 1 CTRL   <tibble [1,245 × 14]> 0.0595 0.113  0.0701 0.135   0.829  0.930
# 2 KI36   <tibble [1,618 × 14]> 0.0626 0.0812 0.0693 0.0833  0.979  0.966
# 3 KIA2   <tibble [1,895 × 14]> 0.0310 0.0324 0.0261 0.0231  0.928  0.917


## plot MMB on UMAP
# OXPHOS complex
lapply(c("CI", "CIII", "CIV", "CV", "tRNA", "rRNA", "Total"), function(cat){
  
  df <- scmtMMB_HEK %>% dplyr::filter(symbol == cat & barcode != "pseudobulk") %>% 
    mutate(log10_MpMB = log10(mutation_per_MB), log2_MSSw = log2(MSS_weighted)) %>% 
    mutate(log10_MpMB = ifelse(is.infinite(log10_MpMB), NA, log10_MpMB),
           log2_MSSw = ifelse(is.infinite(log2_MSSw), NA, log2_MSSw))
  
  p <- ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, color = log10_MpMB)) +
    geom_point(size = 0.01) + 
    scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra"), 
                          limit = quantile(df$log10_MpMB, c(0.02,0.98), na.rm =TRUE)) +
    theme(aspect.ratio=1) + theme_classic() + ggtitle(paste0("MpMB_", cat)) 
  p.nl <- p + theme_void() + Seurat::NoLegend() 
  p.leg <- ggpubr::get_legend(p) %>% ggpubr::as_ggplot()
  ggsave(plot = p.nl, paste0("../plot/5_MpMB_UMAP/MpMB_", cat,".pdf"), width = 2, height = 2)
  ggsave(plot = p.leg, paste0("../plot/5_MpMB_UMAP/MpMB_", cat,".leg.pdf"), width = 3, height = 2)
  
  p <- ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, color = log2_MSSw)) +
    geom_point(size = 0.01) + scale_colour_viridis_c() +
    theme(aspect.ratio=1) + theme_classic() + ggtitle(paste0("MSS_weighted_", cat))
  p.nl <- p + theme_void() + Seurat::NoLegend() 
  p.leg <- ggpubr::get_legend(p) %>% ggpubr::as_ggplot()
  ggsave(plot = p.nl, paste0("../plot/5_MSSw_UMAP/MSSw_", cat,".pdf"), width = 2, height = 2)
  ggsave(plot = p.leg, paste0("../plot/5_MSSw_UMAP/MSSw_", cat,".leg.pdf"), width = 3, height = 2)
  
})

# gene
lapply(unique(scmtMMB_HEK$symbol), function(cat){
  
  df <- scmtMMB_HEK %>% dplyr::filter(symbol == cat & barcode != "pseudobulk") %>% 
    mutate(log10_MpMB = log10(mutation_per_MB), log2_MSSw = log2(MSS_weighted)) %>% 
    mutate(log10_MpMB = ifelse(is.infinite(log10_MpMB), NA, log10_MpMB),
           log2_MSSw = ifelse(is.infinite(log2_MSSw), NA, log2_MSSw))
  
  p <- ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, color = log10_MpMB)) +
    geom_point(size = 0.01) + 
    scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra"), 
                          limit = quantile(df$log10_MpMB, c(0.02,0.98), na.rm =TRUE)) +
    theme(aspect.ratio=1) + theme_classic() + ggtitle(paste0("MpMB_", cat)) 
  p.nl <- p + theme_void() + Seurat::NoLegend() 
  p.leg <- ggpubr::get_legend(p) %>% ggpubr::as_ggplot()
  ggsave(plot = p.nl, paste0("../plot/5a_MpMB_UMAP/MpMB_", cat,".pdf"), width = 2, height = 2)
  ggsave(plot = p.leg, paste0("../plot/5a_MpMB_UMAP/MpMB_", cat,".leg.pdf"), width = 3, height = 2)
  
  p <- ggplot(df, aes(x = atacUMAP_1, y = atacUMAP_2, color = log2_MSSw)) +
    geom_point(size = 0.01) + scale_colour_viridis_c() +
    theme(aspect.ratio=1) + theme_classic() + ggtitle(paste0("MSS_weighted_", cat))
  p.nl <- p + theme_void() + Seurat::NoLegend() 
  p.leg <- ggpubr::get_legend(p) %>% ggpubr::as_ggplot()
  ggsave(plot = p.nl, paste0("../plot/5a_MSSw_UMAP/MSSw_", cat,".pdf"), width = 2, height = 2)
  ggsave(plot = p.leg, paste0("../plot/5a_MSSw_UMAP/MSSw_", cat,".leg.pdf"), width = 3, height = 2)
  
})



p7 <- scmtMMB_HEK %>% filter(symbol == "Total") %>% filter(complete.cases(.)) %>% 
  ggplot(aes(x = log10(mtDNA_depth), y = mutation_per_MB + 1)) +
  geom_point(aes(col = sample)) + facet_wrap(~seurat_clusters) + 
  scale_color_manual(values = color_vec) + 
  theme_bw() + ylab("Mutation per MB (Total)") + theme(legend.position="none") +
  geom_smooth(method = "lm")
ggsave(plot = p7, paste0("../plot/57_CorDepth_MpMB_cluster.pdf"), width=6, height=6)

p7s <- ggscatter(
  scmtMMB_HEK_tl, x = "mtDNA_depth", y = "mutation_per_MB",
  color = "seurat_clusters", 
  add = "reg.line"
) + facet_wrap(~seurat_clusters) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
    label.y = 2750)
ggsave(plot = p7s, paste0("../plot/57s_CorDepth_MpMB_cluster.pdf"), width=6, height=6)

# correlation between complex
scmtMMB_HEK %>% dplyr::select(barcode, sample, symbol, mutation_per_MB) %>% 
  dplyr::filter(grepl("CI|CV|RNA", symbol)) %>% 
  mutate(mutation_per_MB = log10(mutation_per_MB + 1)) %>% 
  tidyr::pivot_wider(names_from = symbol, values_from = mutation_per_MB) %>% 
  dplyr::filter(sample == "KIA2") %>% dplyr::select(-sample, -barcode) %>% 
  dplyr::select(CI, CIII, CIV, CV, tRNA, rRNA) %>% 
  psych::pairs.panels(., digits = 3, pch = 21, bg="#000090", smooth = FALSE, ellipses = FALSE)


scmtMMB_HEK %>% dplyr::select(barcode, sample, symbol, mutation_per_MB) %>% 
  dplyr::filter(grepl("CI|CV|RNA", symbol)) %>% 
  mutate(mutation_per_MB = log10(mutation_per_MB + 1)) %>% 
  tidyr::pivot_wider(names_from = symbol, values_from = mutation_per_MB) %>% 
  dplyr::filter(sample == "KI36") %>% dplyr::select(-sample, -barcode) %>% 
  dplyr::select(CI, CIII, CIV, CV, tRNA, rRNA) %>% 
  psych::pairs.panels(., digits = 3, pch = 21, bg="#960096", smooth = FALSE, ellipses = FALSE)

scmtMMB_HEK %>% dplyr::select(barcode, sample, symbol, mutation_per_MB) %>% 
  dplyr::filter(grepl("CI|CV|RNA", symbol)) %>% 
  mutate(mutation_per_MB = log10(mutation_per_MB + 1)) %>% 
  tidyr::pivot_wider(names_from = symbol, values_from = mutation_per_MB) %>% 
  dplyr::filter(sample == "CTRL") %>% dplyr::select(-c(sample, barcode)) %>%
  select(CI, CIII, CIV, CV, tRNA, rRNA) %>% 
  psych::pairs.panels(., digits = 3, pch = 21, bg="#2780FF", smooth = FALSE, ellipses = FALSE)

# correlation between gene

df.list <- scmtMMB_HEK %>% select(barcode, sample, symbol, mutation_per_MB) %>% 
  dplyr::filter(!grepl("CI|CV|RNA|Total", symbol)) %>%
  pivot_wider(names_from = symbol, values_from = mutation_per_MB) %>% 
  select(-barcode) %>% 
  group_split(sample, .keep = FALSE)

names(df.list) <- c("CTRL", "KI36", "KIA2")
df.list <- lapply(df.list, function(M){as.matrix(M) %>% replace(!is.finite(.), 0)})
cor.list <- lapply(df.list, function(M){cor(M, use = "complete.obs")})  
cor.list <- lapply(cor.list, function(M){
  M <- M %>% replace(!is.finite(.), 0)
  colnames(M) <- gsub(pattern = "MT_", replacement = "", colnames(M))
  rownames(M) <- gsub(pattern = "MT_", replacement = "", rownames(M))
  M
})

corrplot::corrplot(cor.list$KIA2, order = 'hclust', type = 'upper', tl.col="black",diag = FALSE, insig = 'label_sig')
corrplot::corrplot(cor.list$KI36, order = 'hclust', type = 'upper', tl.col="black",diag = FALSE, insig = 'label_sig')
corrplot::corrplot(cor.list$CTRL, order = 'hclust', type = 'upper', tl.col="black",diag = FALSE, insig = 'label_sig')

scmtMMB_HEK %>% 
  dplyr::filter(grepl("Total", symbol)) %>% 
  group_by(sample) %>% 
  summarize(across(c(mutation_per_MB, MSS_weighted), mean))

scmtMMB_HEK %>% 
  dplyr::filter(grepl("CI|CV|RNA|Total", symbol)) %>% 
  group_by(sample, symbol) %>% 
  summarize(across(c(mutation_per_MB, MSS_weighted), median)) %>% 
  print(21)