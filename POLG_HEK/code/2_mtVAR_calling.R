setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code")
source("../../global_func/variant_calling.R")
source("../../global_func/mtGene_color.R")
library(Matrix)
library(Seurat)
library(Signac)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(tibble)
library(ggpubr)
library(rstatix)
'%ni%' <- Negate('%in%')

# Load the seurat object
HEK <- readRDS("../POLG_HEK_large_data_files/output/1_POLG_HEK_preprocesssed_seurat.rds")
color_vec <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")

# extract cell barcodes for each line
CTRL_CB <- WhichCells(HEK, idents = "CTRL")
KI36_CB <- WhichCells(HEK, idents = "KI36")
KIA2_CB <- WhichCells(HEK, idents = "KIA2")
CTRL <- subset(HEK, cells = CTRL_CB)
KI36 <- subset(HEK, cells = KI36_CB)
KIA2 <- subset(HEK, cells = KIA2_CB)

# Load the mitochondrial genome data
mito.SE <- readRDS("../POLG_HEK_large_data_files/input/mgatk/mgatk.rds")[,Cells(HEK)]
cov_per_pos_cell <- assays(mito.SE)[['coverage']]
mitoSE.CTRL <- mito.SE[,CTRL_CB]
mitoSE.KI36 <- mito.SE[,KI36_CB]
mitoSE.KIA2 <- mito.SE[,KIA2_CB]

# save RDS file
saveRDS(cov_per_pos_cell[,CTRL_CB], "../POLG_HEK_large_data_files/output/2_mtDNA_pos_coverage.CTRL.rds")
saveRDS(cov_per_pos_cell[,KI36_CB], "../POLG_HEK_large_data_files/output/2_mtDNA_pos_coverage.KI36.rds")
saveRDS(cov_per_pos_cell[,KIA2_CB], "../POLG_HEK_large_data_files/output/2_mtDNA_pos_coverage.KIA2.rds")
saveRDS(mitoSE.CTRL, "../POLG_HEK_large_data_files/output/2_mgatk.CTRL.rds")
saveRDS(mitoSE.KI36, "../POLG_HEK_large_data_files/output/2_mgatk.KI36.rds")
saveRDS(mitoSE.KIA2, "../POLG_HEK_large_data_files/output/2_mgatk.KIA2.rds")

## Visualize position-wise mitochondral genome coverage
pull_coverage <- function(SE, cells, resolution = 5){
  zoo::rollmean(rowMeans(assays(SE)[['coverage']]), resolution)
}

cov_df <- data.frame(
  pos = zoo::rollmean(1:16569, 5),
  CTRL = pull_coverage(mitoSE.CTRL),
  KI36 = pull_coverage(mitoSE.KI36),
  KIA2 = pull_coverage(mitoSE.KIA2))

p1 <- cov_df %>% reshape2::melt(id.vars = "pos") %>% # dplyr::filter(value > 5) %>%
  ggplot() +
  geom_line(aes(x = pos, y = value, color = variable)) +
  scale_y_continuous(limits = c(-3, 75), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,16596, by = 5000), expand = c(0,0)) + 
  scale_color_manual(values = color_vec) +
  labs(x = "Position on mtDNA chromosome", y = "Roll mean coverage") + geom_rect(data = mitochondrial_data, 
                                                                                 aes(xmin = Starting, xmax = Ending, ymin = -3, ymax = 0, fill = color), alpha = 1) +
  scale_fill_manual(values = mtGene_colors) +
  # coord_polar() + 
  guides(fill=guide_legend(ncol=2))+
  theme_classic()

p1.nl <- p1 + Seurat::NoLegend()
leg <- ggpubr::get_legend(p1) %>% ggpubr::as_ggplot()
ggsave("../plot/21_mtDNA_cov.pdf", p1.nl, width = 5, height = 3)
ggsave("../plot/21_mtDNA_cov_leg.pdf", leg, width = 5, height = 5)

## Mitochondrial variant calling
mmat.CTRL <- call_mutations_mgatk(mitoSE.CTRL)
mmat.KI36 <- call_mutations_mgatk(mitoSE.KI36)
mmat.KIA2 <- call_mutations_mgatk(mitoSE.KIA2)

VarPlot <- function(mmat, SampleID){
  counts <- as.data.frame(rowData(mmat)) %>% subset(n_cells_conf_detected >= 1 & 
                                                      strand_correlation > 0.65 & vmr > 0.01) %>% nrow()
  plt <- VariantPlot(as.data.frame(rowData(mmat))) + 
    labs(title = SampleID, subtitle = paste0("Number of somatic variants: ", counts)) + 
    scale_y_log10(limits = c(-5,1), n.breaks = 7) + scale_x_continuous(limits = c(-0.6,1), n.breaks = 4)
  ggsave(plot = plt, paste0("../plot/22_VMR_",SampleID,".pdf"), width = 3, height = 3)
  return(plt)
}

p2.1 <- VarPlot(mmat.CTRL, "CTRL")
p2.2 <- VarPlot(mmat.KI36, "KI36")
p2.3 <- VarPlot(mmat.KIA2, "KIA2")
cowplot::plot_grid(p2.1, p2.2, p2.3, nrow = 1)

### mtDNA somatic variant
sumVar <- function(mmat, n_cells = 1){
  var <- as.data.frame(rowData(mmat))
  var <- var %>% subset(n_cells_conf_detected >= n_cells &
                          variant %ni% c("301A>C", "302A>C", "310T>C", "316G>C")) %>% 
    mutate(type = ifelse(strand_correlation > 0.65,
                         ifelse(vmr > 0.01, "high-confident", "germline"), "error-prone"))
  table(var$type)
}

sumVar(mmat.CTRL)
sumVar(mmat.KI36)
sumVar(mmat.KIA2)

mtMut <- function(mmat, n_cells = 1, vmr.thres = 0.01){
  var <- as.data.frame(rowData(mmat))
  var <- var %>% subset(n_cells_conf_detected >= n_cells &
                          strand_correlation > 0.65 &
                          vmr >= vmr.thres &
                          variant %ni% c("301A>C", "302A>C", "310T>C", "316G>C")) %>% pull(variant)
  mmat <- mmat[var,]
  return(mmat)
}

# analyze three cell lines separately and visualize selected variants
mmatMut.CTRL <- mtMut(mmat.CTRL)
rowData(mmatMut.CTRL)$condition <- "CTRL"
mmatMut.KI36 <- mtMut(mmat.KI36)
rowData(mmatMut.KI36)$condition <- "KI36"
mmatMut.KIA2 <- mtMut(mmat.KIA2)
rowData(mmatMut.KIA2)$condition <- "KIA2"


saveRDS(mmatMut.CTRL, "../POLG_HEK_large_data_files/output/2_mmtx.CTRL.rds")
saveRDS(mmatMut.KI36, "../POLG_HEK_large_data_files/output/2_mmtx.KI36.rds")
saveRDS(mmatMut.KIA2, "../POLG_HEK_large_data_files/output/2_mmtx.KIA2.rds")
mmatMut.CTRL <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.CTRL.rds")
mmatMut.KI36 <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KI36.rds")
mmatMut.KIA2 <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KIA2.rds")

df1 <- data.frame(
  condition = c(rep(x = "CTRL", ncol(mmatMut.CTRL)), rep(x = "KI36", ncol(mmatMut.KI36)), rep(x = "KIA2", ncol(mmatMut.KIA2))),
  counts = c(colSums(assay(mmatMut.CTRL)>0), colSums(assay(mmatMut.KI36)>0), colSums(assay(mmatMut.KIA2)>0)),
  mtDNA_depth = c(mmatMut.CTRL$depth, mmatMut.KI36$depth, mmatMut.KIA2$depth)) %>% 
  mutate(counts_norm = counts/mtDNA_depth) 


p <- df1 %>% 
  ggplot(aes(x= condition, y = counts)) +
  geom_violin(scale = "width", aes(fill = condition)) + geom_boxplot(width=0.1, outlier.size=0) + 
  scale_fill_manual(values = color_vec) + ylab("mtDNA mutation counts") +
  theme_classic() + NoLegend()

ggsave(plot = p, "../plot/2_mtVAR_counts.pdf", width = 2.5, height = 2.5)

ggplot(data = df1, aes(x=mtDNA_depth, y = counts, color = condition)) +
  geom_point(size=0.01) + scale_color_manual(values = color_vec) + geom_smooth(method = "lm") +
  facet_wrap(~condition) +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) + 
  theme_classic() + NoLegend()

pn <- df1 %>% 
  ggplot(aes(x= condition, y = counts/mtDNA_depth)) +
  geom_violin(scale = "width", aes(fill = condition)) + geom_boxplot(width=0.1, outlier.size=0) + 
  scale_fill_manual(values = color_vec) + ylab("Normalized mtDNA mutation counts") +
  theme_classic() + NoLegend()

ggsave(plot = pn, "../plot/2_normalized_mtVAR_counts.pdf", width = 2.5, height = 2.5)

summarizer <- function(data, numeric_cols = NULL, ...) {
  data %>%
    group_by(...) %>%
    summarise(across({{numeric_cols}}, list(
      min = ~min(.x),
      q05 = ~quantile(.x, 0.05, na.rm = TRUE),
      median = ~median(.x, na.rm = TRUE),
      q95 = ~quantile(.x, 0.95, na.rm = TRUE),
      max = ~max(.x)
    ), .names = "{col}_{fn}"))
}

summarizer(df1, c(counts, counts_norm), condition) %>% tidyr::pivot_longer(cols = -condition, names_to = "name", values_to = "value") %>% 
  mutate(
    feature = sub("^(.*)_[^_]+$", "\\1", name),
    stat = sub("^.*_([^_]+)$", "\\1", name) 
  ) %>% 
  select(condition, feature, stat, value) %>% 
  tidyr::pivot_wider(names_from = condition, values_from = value)

#    feature     stat      CTRL   KI36    KIA2
#    <chr>       <chr>    <dbl>  <dbl>   <dbl>
#  1 counts      min      2      51     110   
#  2 counts      q05      4     103.    254.  
#  3 counts      median   9     180     472.  
#  4 counts      q95     20     288     721   
#  5 counts      max    218     592    1123   
#  6 counts_norm min      0.146   2.38    5.12
#  7 counts_norm q05      0.353   3.91    7.41
#  8 counts_norm median   0.769   6.38   10.9 
#  9 counts_norm q95      1.98   10.5    17.1 
# 10 counts_norm max      5.34   15      25.5 


## mtDNA somatic variant visualization 

CTRL[["allele"]] <- CreateAssayObject(assays(mmatMut.CTRL)$allele_frequency)
DefaultAssay(CTRL) <- "allele"
px1 <-  DoHeatmap(CTRL, slot = "data", disp.max = 1, group.colors = color_vec,
                  features = as.data.frame(rowData(mmatMut.CTRL)) %>% 
                    arrange(desc(mean), desc(n_cells_conf_detected)) %>%
                    slice_head(n = 10) %>%  pull(variant)) + scale_fill_viridis_c() + NoLegend() 

KI36[["allele"]] <- CreateAssayObject(assays(mmatMut.KI36)$allele_frequency)
DefaultAssay(KI36) <- "allele"
px2 <- DoHeatmap(KI36, slot = "data", disp.max = 1, group.colors = color_vec,
                 features = as.data.frame(rowData(mmatMut.KI36)) %>% 
                   arrange(desc(mean), desc(n_cells_conf_detected)) %>%
                   slice_head(n = 10) %>%  pull(variant)) + scale_fill_viridis_c() + NoLegend() 
# labs(title = "Top 10 high mean heteroplamy variants")

KIA2[["allele"]] <- CreateAssayObject(assays(mmatMut.KIA2)$allele_frequency)
DefaultAssay(KIA2) <- "allele"
px3 <- DoHeatmap(KIA2, slot = "data", disp.max = 1, group.colors = color_vec,
                 features = as.data.frame(rowData(mmatMut.KIA2)) %>% 
                   arrange(desc(mean), desc(n_cells_conf_detected)) %>%
                   slice_head(n = 10) %>%  pull(variant)) + scale_fill_viridis_c() + NoLegend() 
#labs(title = "Top 10 high mean heteroplamy variants")

cowplot::plot_grid(px1, px2, px3, nrow = 1)

# **Unique and share somatic mutations**
Venn <- BioVenn::draw.venn(
  rownames(mmatMut.CTRL), rownames(mmatMut.KI36), rownames(mmatMut.KIA2),
  title = NULL, subtitle = NULL,
  xtitle = NULL, ytitle = NULL, ztitle = NULL, 
  nr_c = "black", x_c = "#2780FF", y_c = "#960096", z_c = "#000090",
  bg_c = "white", output = "pdf", filename = "../plot/23_mtVar_venn.pdf", width = 800, height =800
)

# [1] "x total: 455", "y total: 2483", "z total: 4651"
# [1] "x only: 67", "y only: 1576", "z only: 3569"
# [1] "x-y total overlap: 161", "x-z total overlap: 336", "y-z total overlap: 855"
# [1] "x-y only overlap: 52", "x-z only overlap: 227", "y-z only overlap: 746", "x-y-z overlap: 109"

# **Heteroplasmy stratified by unique and share somatic mutations**
mtVar.bulk.df <- do.call(rbind, lapply(list(rowData(mmatMut.CTRL), rowData(mmatMut.KI36), rowData(mmatMut.KIA2)), as_tibble))
mtVar_shared <- Venn$xyz; mtVar_WT_only <- Venn$x_only; mtVar_KI36_only <- Venn$y_only; mtVar_KIA2_only <- Venn$z_only
# mtVar_parental <- Venn$xyz; mtVar_WT_sahred <- c(Venn$xy_only, Venn$xz_only)
mtVar_parental <- c(Venn$xy_only, Venn$xz_only, Venn$xyz)
mtVar_POLG_sahred <- Venn$yz_only; mtVar_POLG_only <- c(mtVar_KI36_only, mtVar_POLG_sahred, mtVar_KIA2_only)

mtVar.bulk.df <- mtVar.bulk.df %>% 
  mutate(occurrence = case_when(variant %in% mtVar_parental ~ "Parental",
                                variant %in% mtVar_POLG_sahred ~ "POLG-shared",
                                TRUE ~ "de-novo")) %>% 
  mutate(occurrence = factor(occurrence, levels = c("Parental", "POLG-shared", "de-novo")))

write.csv(mtVar.bulk.df, file = "../output/2_mtVar.bulk.meta.csv", 
          quote = FALSE, row.names = TRUE)


p4 <- mtVar.bulk.df %>% ggplot(aes(x = condition, y = mean)) +
  geom_violin(aes(fill = condition), width = 0.7) + geom_boxplot(width=0.1) + theme_classic() + 
  scale_fill_manual(values = color_vec) + scale_y_log10(limits = c(5e-5,1), breaks = 10^(-5:0)) + Seurat::NoLegend()
ggsave(plot = p4 , paste0("../plot/24_heteroplasmy_bulk.pdf"), width = 3, height = 3) 


stat.test <- mtVar.bulk.df %>% 
  wilcox_test(mean ~ condition) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
write.csv(stat.test, file = "../output/2_het_wilcox.test.csv", 
          quote = FALSE, row.names = TRUE)

lapply(c("CTRL"), function(cdt){
  colors =  c("Parental" = "#E8C748", "POLG-shared" = "#583690","Private" = as.character(color_vec[cdt])) #"WT-shared" = "#1066FF"
  p5 <- mtVar.bulk.df %>% filter(condition == cdt) %>% ggplot(aes(x = occurrence, y = mean)) +
    geom_violin(aes(fill = occurrence), width = 0.7, scale = "width") + geom_boxplot(width=0.1) + theme_classic() + 
    scale_fill_manual(values = colors) + scale_y_log10(limits = c(5e-5,1), breaks = 10^(-5:0)) +
    ggtitle(cdt) + Seurat::NoLegend()
  print(p5)
  ggsave(plot = p5, paste0("../plot/25_heteroplasmy_bulk_occur_",cdt,".pdf"), width = 2.3, height = 2)
})

lapply(c("CTRL", "KI36","KIA2"), function(cdt){
  colors =  c("Parental" = "#E8C748", "POLG-shared" = "#583690","de-novo" = as.character(color_vec[cdt])) #"WT-shared" = "#1066FF"
  p5 <- mtVar.bulk.df %>% filter(condition == cdt) %>% ggplot(aes(x = occurrence, y = mean)) +
    geom_violin(aes(fill = occurrence), width = 0.7, scale = "width") + geom_boxplot(width=0.1) + theme_classic() + 
    scale_fill_manual(values = colors) + scale_y_log10(limits = c(5e-5,1), breaks = 10^(-5:0)) +
    ggtitle(cdt) + Seurat::NoLegend()
  print(p5)
  ggsave(plot = p5, paste0("../plot/25_heteroplasmy_bulk_occur_",cdt,".pdf"), width = 3, height = 2)
})

stat.test.occu <- mtVar.bulk.df %>% select(condition, occurrence, mean) %>% 
  group_by(condition) %>% 
  wilcox_test(mean ~ occurrence) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test.occu
write.csv(stat.test.occu, file = "../output/2_het_wilcox.test.occu.csv", 
          quote = FALSE, row.names = TRUE)

lapply(c("Parental", "de-novo","POLG-shared"), function(occu){
  p6 <- mtVar.bulk.df %>% filter(occurrence == occu) %>% ggplot(aes(x = condition, y = mean)) +
    geom_violin(aes(fill = condition), width = 0.7, scale = "width") + geom_boxplot(width=0.1) + theme_classic() + 
    scale_fill_manual(values = color_vec) + scale_y_log10(limits = c(5e-5,1), breaks = 10^(-5:0)) +
    ggtitle(occu)+ Seurat::NoLegend()
  print(p6)
  ggsave(plot = p6, paste0("../plot/26_heteroplasmy_bulk_cdt_",occu,".pdf"), width = 3, height = 2)
})

stat.test.cdt <- mtVar.bulk.df %>% select(condition, occurrence, mean) %>% 
  group_by(occurrence) %>% 
  wilcox_test(mean ~ condition) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test.cdt
write.csv(stat.test.cdt, file = "../output/2_het_wilcox.test.cdt.csv", 
          quote = FALSE, row.names = TRUE)

mtVar.bulk.df %>% filter(occurrence == "Parental") %>% ggplot(aes(x = condition, y = mean, group = variant)) + 
  geom_line() + theme_classic() + 
  scale_fill_manual(values = color_vec) + scale_y_log10(limits = c(5e-5,1), breaks = 10^(-5:0)) 


p7 <- mtVar.bulk.df %>% filter(occurrence != "Parental") %>% ggplot(aes(x = condition, y = mean)) +
  geom_violin(aes(fill = condition), width = 0.7, scale = "width") + geom_boxplot(width=0.1) + theme_classic() + 
  scale_fill_manual(values = color_vec) + scale_y_log10(limits = c(5e-5,1), breaks = 10^(-5:0)) +
  ggtitle("all de novo")+ Seurat::NoLegend()
print(p7)
ggsave(plot = p7, paste0("../plot/27_heteroplasmy_bulk_cdt_all_de_novo.pdf"), width = 3, height = 2)

stat.test.all.denovo <- mtVar.bulk.df %>% select(condition, occurrence, mean) %>% 
  mutate(occurrence = case_when(occurrence == "Parental" ~ "Parental", TRUE ~ "de novo")) %>% 
  group_by(occurrence) %>% 
  wilcox_test(mean ~ condition) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test.all.denovo
write.csv(stat.test.all.denovo, file = "../output/2_het_wilcox.test.all.denovo.csv", 
          quote = FALSE, row.names = TRUE)






