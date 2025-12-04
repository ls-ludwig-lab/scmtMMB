setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
library(Matrix)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(tibble)
library(ggpubr)
library(rstatix)
'%ni%' <- Negate('%in%')


mmatMut.CTRL <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.CTRL.rds")
mmatMut.KI36 <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KI36.rds")
mmatMut.KIA2 <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KIA2.rds")

# Load the seurat object
color_vec <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")
df1 <- data.frame(
  condition = c(rep(x = "CTRL", ncol(mmatMut.CTRL)), rep(x = "KI36", ncol(mmatMut.KI36)), rep(x = "KIA2", ncol(mmatMut.KIA2))),
  counts = c(colSums(assay(mmatMut.CTRL)>0), colSums(assay(mmatMut.KI36)>0), colSums(assay(mmatMut.KIA2)>0)),
  mtDNA_depth = c(mmatMut.CTRL$depth, mmatMut.KI36$depth, mmatMut.KIA2$depth)) %>% 
  mutate(counts_norm = counts/mtDNA_depth) 


p <- df1 %>% 
  ggplot(aes(x= condition, y = counts)) +
  geom_violin(scale = "width", aes(fill = condition)) + geom_boxplot(width=0.1, outlier.size=0) + 
  scale_fill_manual(values = color_vec) + ylab("mtDNA mutation counts") +
  theme_classic() + theme(legend.position = "none") + 
  scale_y_log10(breaks = c(10, 100, 1000, 3000), 
                labels = c(10, 100, 1000, 3000))

ggsave(plot = p, "../plot/Fig2a_2_mtVAR_counts.pdf", width = 2.5, height = 2.5) # Fig2a

ggplot(data = df1, aes(x=mtDNA_depth, y = counts, color = condition)) +
  geom_point(size=0.01) + scale_color_manual(values = color_vec) + geom_smooth(method = "lm") +
  facet_wrap(~condition) +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) + 
  theme_classic() + theme(legend.position = "none")

pn <- df1 %>% 
  ggplot(aes(x= condition, y = counts/mtDNA_depth)) +
  geom_violin(scale = "width", aes(fill = condition)) + geom_boxplot(width=0.1, outlier.size=0) + 
  scale_fill_manual(values = color_vec) + ylab("Normalized mtDNA mutation counts") +
  theme_classic() + theme(legend.position = "none")

ggsave(plot = pn, "../plot/Fig2b_2_normalized_mtVAR_counts.pdf", width = 2.5, height = 2.5) # Fig2b

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

# feature     stat      CTRL    KI36    KIA2
# <chr>       <chr>    <dbl>   <dbl>   <dbl>
# 1 counts      min      8      193     229   
# 2 counts      q05     18      488     844   
# 3 counts      median  35      904    1487   
# 4 counts      q95     62     1468.   2106.  
# 5 counts      max    488     2328    3386   
# 6 counts_norm min      0.184    3.64    3.36
# 7 counts_norm q05      0.327    4.68    4.95
# 8 counts_norm median   0.750    7.01    6.97
# 9 counts_norm q95      2.23    12.7    12.5 
# 10 counts_norm max      4.76    20.5    23.7 


# **Unique and share somatic mutations**
Venn <- BioVenn::draw.venn(
  rownames(mmatMut.CTRL), rownames(mmatMut.KI36), rownames(mmatMut.KIA2),
  title = NULL, subtitle = NULL,
  xtitle = NULL, ytitle = NULL, ztitle = NULL, 
  nr_c = "black", x_c = "#2780FF", y_c = "#960096", z_c = "#000090",
  bg_c = "white", output = "pdf", filename = "../plot/Fig2c_23_mtVar_venn.pdf", width = 800, height =800
)

# [1] "x total: 620" [1] "y total: 9655" [1] "z total: 11407"
# [1] "x only: 27" [1] "y only: 4125" [1] "z only: 5698"
# [1] "x-y total overlap: 407" [1] "x-z total overlap: 586" [1] "y-z total overlap: 5523"
# [1] "x-y only overlap: 7" [1] "x-z only overlap: 186" [1] "y-z only overlap: 5123" [1] "x-y-z overlap: 400"

# **Heteroplasmy stratified by unique and share somatic mutations**
mtVar.bulk.df <- do.call(rbind, lapply(list(rowData(mmatMut.CTRL), rowData(mmatMut.KI36), rowData(mmatMut.KIA2)), as_tibble))
mtVar_shared <- Venn$xyz; mtVar_WT_only <- Venn$x_only; mtVar_KI36_only <- Venn$y_only; mtVar_KIA2_only <- Venn$z_only
# mtVar_parental <- Venn$xyz; mtVar_WT_sahred <- c(Venn$xy_only, Venn$xz_only)
mtVar_parental <- c(Venn$xy_only, Venn$xz_only, Venn$xyz)
mtVar_POLG_sahred <- Venn$yz_only; mtVar_POLG_only <- c(mtVar_KI36_only, mtVar_POLG_sahred, mtVar_KIA2_only)

mtVar.bulk.df <- mtVar.bulk.df %>% 
  mutate(occurrence = case_when(variant %in% mtVar_parental ~ "Parental",
                                variant %in% mtVar_POLG_sahred ~ "POLG-shared",
                                TRUE ~ "Line-specific")) %>% 
  mutate(occurrence = factor(occurrence, levels = c("Parental", "POLG-shared", "Line-specific")))

write.csv(mtVar.bulk.df, file = "../output/2_mtVar.bulk.meta.csv", 
          quote = FALSE, row.names = TRUE)


p4 <- mtVar.bulk.df %>% ggplot(aes(x = condition, y = mean)) +
  geom_violin(aes(fill = condition), width = 0.7) + geom_boxplot(width=0.1) + theme_classic() + 
  scale_fill_manual(values = color_vec) + scale_y_log10(limits = c(1e-5,1), breaks = 10^(-5:0)) + Seurat::NoLegend()
ggsave(plot = p4 , paste0("../plot/24_heteroplasmy_bulk.pdf"), width = 3, height = 3) 


stat.test <- mtVar.bulk.df %>% 
  wilcox_test(mean ~ condition) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
stat.test
write.csv(stat.test, file = "../output/2_het_wilcox.test.csv", 
          quote = FALSE, row.names = TRUE)

lapply(c("CTRL"), function(cdt){
  colors =  c("Parental" = "#E8C748", "POLG-shared" = "#583690", "Line-specific" = as.character(color_vec[cdt])) #"WT-shared" = "#1066FF"
  p5 <- mtVar.bulk.df %>% filter(condition == cdt) %>% ggplot(aes(x = occurrence, y = mean)) +
    geom_violin(aes(fill = occurrence), width = 0.7, scale = "width") + geom_boxplot(width=0.1) + theme_classic() + 
    scale_fill_manual(values = colors) + scale_y_log10(limits = c(1e-5,1), breaks = 10^(-5:0)) +
    ggtitle(cdt) + Seurat::NoLegend()
  print(p5)
  ggsave(plot = p5, paste0("../plot/Fig2g_25_heteroplasmy_bulk_occur_",cdt,".pdf"), width = 2, height = 2)
})


lapply(c("KI36","KIA2"), function(cdt){
  colors =  c("Parental" = "#E8C748", "POLG-shared" = "#583690", "Line-specific" = as.character(color_vec[cdt])) #"WT-shared" = "#1066FF"
  p5 <- mtVar.bulk.df %>% filter(condition == cdt) %>% ggplot(aes(x = occurrence, y = mean)) +
    geom_violin(aes(fill = occurrence), width = 0.7, scale = "width") + geom_boxplot(width=0.1) + theme_classic() + 
    scale_fill_manual(values = colors) + scale_y_log10(limits = c(1e-5,1), breaks = 10^(-5:0)) +
    ggtitle(cdt) + Seurat::NoLegend()
  print(p5)
  ggsave(plot = p5, paste0("../plot/Fig2g_25_heteroplasmy_bulk_occur_",cdt,".pdf"), width = 3, height = 2)
})

stat.test.occu <- mtVar.bulk.df %>% select(condition, occurrence, mean) %>% 
  group_by(condition) %>% 
  wilcox_test(mean ~ occurrence) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
stat.test.occu
# condition .y.   group1      group2           n1    n2 statistic         p     p.adj p.adj.signif
# <chr>     <chr> <chr>       <chr>         <int> <int>     <dbl>     <dbl>     <dbl> <chr>       
#   1 CTRL      mean  Parental    Line-specific   593    27    12298  2.42e-  6 2.42e-  6 ****        
#   2 KI36      mean  Parental    POLG-shared     407  5123  1390815  2.75e- 29 3.21e- 29 ****        
#   3 KI36      mean  Parental    Line-specific   407  4125  1288696. 3.48e- 71 4.87e- 71 ****        
#   4 KI36      mean  POLG-shared Line-specific  5123  4125 12941331  2.64e- 77 4.62e- 77 ****        
#   5 KIA2      mean  Parental    POLG-shared     586  5123  2343685  4.14e-110 9.66e-110 ****        
#   6 KIA2      mean  Parental    Line-specific   586  5698  2937316. 6.95e-202 2.43e-201 ****        
#   7 KIA2      mean  POLG-shared Line-specific  5123  5698 20020800. 3.82e-245 2.67e-244 ****    
write.csv(stat.test.occu, file = "../output/2_het_wilcox.test.occu.csv", 
          quote = FALSE, row.names = TRUE)

lapply(c("Parental", "POLG-shared"), function(occu){
  p6 <- mtVar.bulk.df %>% filter(occurrence == occu) %>% ggplot(aes(x = condition, y = mean)) +
    geom_violin(aes(fill = condition), width = 0.7, scale = "width") + geom_boxplot(width=0.1) + theme_classic() + 
    scale_fill_manual(values = color_vec) + scale_y_log10(limits = c(1e-5,1), breaks = 10^(-5:0)) +
    ggtitle(occu)+ Seurat::NoLegend()
  print(p6)
  ggsave(plot = p6, paste0("../plot/Fig2f_6_heteroplasmy_bulk_cdt_",occu,".pdf"), width = 2, height = 2)
})

stat.test.cdt <- mtVar.bulk.df %>% select(condition, occurrence, mean) %>% 
  group_by(occurrence) %>% 
  wilcox_test(mean ~ condition) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
stat.test.cdt

# occurrence      .y.   group1 group2     n1    n2 statistic        p    p.adj  p.adj.signif
# <fct>             <chr> <chr>  <chr>  <int> <int>     <dbl>    <dbl>    <dbl> <chr>       
#   1 Parental      mean  CTRL   KI36     593   407    75675  1.14e-23 2.66e-23 ****        
#   2 Parental      mean  CTRL   KIA2     593   586    62234. 3.93e-81 2.75e-80 ****        
#   3 Parental      mean  KI36   KIA2     407   586    72491  6.98e-26 2.44e-25 ****        
#   4 POLG-shared   mean  KI36   KIA2    5123  5123 12053672  9.33e-13 1.31e-12 ****        
#   5 Line-specific mean  CTRL   KI36      27  4125    30945  6.74e- 5 7.86e- 5 ****        
#   6 Line-specific mean  CTRL   KIA2      27  5698    49098. 1   e- 3 1   e- 3 ***         
#   7 Line-specific mean  KI36   KIA2    4125  5698 12916038  4.84e-17 8.47e-17 ****        
  
write.csv(stat.test.cdt, file = "../output/2_het_wilcox.test.cdt.csv", 
          quote = FALSE, row.names = TRUE)


p7 <- mtVar.bulk.df %>% filter(occurrence != "Parental") %>% ggplot(aes(x = condition, y = mean)) +
  geom_violin(aes(fill = condition), width = 0.7, scale = "width") + geom_boxplot(width=0.1) + theme_classic() + 
  scale_fill_manual(values = color_vec) + scale_y_log10(limits = c(1e-5,1), breaks = 10^(-5:0)) +
  ggtitle("Line-specific")+ Seurat::NoLegend()
print(p7)
ggsave(plot = p7, paste0("../plot/Fig2f_27_heteroplasmy_bulk_cdt_Line-specific.pdf"), width = 2, height = 2)

stat.test.all.denovo <- mtVar.bulk.df %>% select(condition, occurrence, mean) %>% 
  mutate(occurrence = case_when(occurrence == "Parental" ~ "Parental", TRUE ~ "Non-parental")) %>% 
  group_by(occurrence) %>% 
  wilcox_test(mean ~ condition) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test.all.denovo
write.csv(stat.test.all.denovo, file = "../output/2_het_wilcox.test.all.denovo.csv", 
          quote = FALSE, row.names = TRUE)






