setwd("~/Ludwig_lab/scmtMMB/PBMC/code/")
source("../../global_func/quantifyMMB.R")
library(Matrix)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(ggplot2)
library(scales)
library(data.table)
library(ggpubr)
library(purrr)
'%ni%' <- Negate('%in%')

## Calculate MMB
color_vec <- c("M80" = "darkorchid4","M60" = "firebrick4", "M35" = "violetred3", "M29" = "firebrick1", "H47" = "blue3", "H05" = "dodgerblue3")

# load metadata
md <- do.call(rbind,lapply(names(color_vec),function(SID){
  df <- fread(paste0("../output/01_metadata_",SID,".tsv"), header = TRUE, sep = "\t")
  mmat <- readRDS(paste0("../PBMC_large_data_files/output/01_mmat_", SID,".rds"))
  if(grepl("H", SID)){
    df$X3243A_G <- 0
  } else{
    df$X3243A_G <- assay(mmat)["3243A>G", df$barcode]
  }
  df
}))

# Calculate MMB exclude 3243AG and the highest VAF variant

scmtMMB_exMax <- do.call(rbind,lapply(names(color_vec),function(SID){
  mmat <- readRDS(paste0("../PBMC_large_data_files/output/01_mmat_", SID,".rds"))
  mmat <- mmat[rownames(mmat) != "3243A>G", ]
  
  # Use apply to modify each column
  mtx <- apply(assay(mmat), 2, function(col) {
    col[which.max(col)] <- 0  # Set max value in each column to 0
    return(col)
  })
  
  # Convert back to matrix if needed
  assay(mmat) <- as(mtx, "dgCMatrix")
  
  rowData(mmat)$mean = rowSums(assays(mmat)$allele_frequency * assays(mmat)$coverage)/rowSums(assays(mmat)$coverage)
  
  cov <- readRDS(paste0("../PBMC_large_data_files/output/01_mtDNA_pos_coverage_", SID,".rds"))
  rbind(
    quantify_MMB(mmat, cov, SID),
    quantify_MMB_psedobulk(mmat, cov, SID)
  )
}))


scmtMMB_exMax %>% mutate(age = as.numeric(gsub("H|M", "", sample))) %>% filter(symbol == "Total") %>% 
  group_by(sample) %>%
  mutate(
    q5 = quantile(mutation_per_MB, 0.05),
    q95 = quantile(mutation_per_MB, 0.95)
  ) %>%
  filter(mutation_per_MB >= q5, mutation_per_MB <= q95) %>%  # Filter out values outside 5th and 95th percentiles
  ungroup() %>% 
  ggplot(aes(x= age, y = mutation_per_MB)) +
  geom_violin(aes(fill = sample), scale = "width") + scale_fill_manual(values = color_vec) + 
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), shape = 23, size = 1, fill = "white") +
  Seurat::NoLegend() + geom_smooth(method = "lm") + theme_bw() + 
  theme(aspect.ratio=1/1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Total scMPM") + xlab("Age (Years)") + Seurat::NoLegend()

scmtMMB_exMax <- left_join(scmtMMB_exMax, md %>% select(predicted.celltype.l1, predicted.celltype.l2, 
                                            refUMAP_1, refUMAP_2, SID, mtDNA_depth, barcode, X3243A_G), 
                     by = c("barcode", "sample" = "SID"))

write.csv(scmtMMB_exMax, "../PBMC_large_data_files/output/2_scmtMMB_exMax.csv", quote = FALSE, row.names = TRUE)

scmtMMB <- read.csv("../PBMC_large_data_files/output/2_scmtMMB.csv", row.names = "X")
scmtMMB$ex <- "all"
scmtMMB_ex <- read.csv("../PBMC_large_data_files/output/2_scmtMMB_exclude_3243.csv", row.names = "X") 
scmtMMB_ex$ex <- "-3243"
scmtMMB_exMax <- read.csv("../PBMC_large_data_files/output/2_scmtMMB_exMax.csv", row.names = "X")
scmtMMB_exMax$ex <- "-3243-max"

color_vec <- c("M80" = "darkorchid4","M60" = "firebrick4", "M35" = "violetred3", 
               "M29" = "firebrick1", "H47" = "blue3", "H05" = "dodgerblue3")

create_violin_plot <- function(data, y_var, y_label, pseudobulk_filter, plot_file) {
  plot <- data %>%
    filter(symbol == "Total") %>%
    mutate(X3243 = factor(ex, levels = c("all", "-3243", "-3243-max"))) %>%
    group_by(sample, ex) %>%
    mutate(
      q5 = quantile(.data[[y_var]], 0.05),
      q95 = quantile(.data[[y_var]], 0.95)
    ) %>%
    filter(.data[[y_var]] >= q5, .data[[y_var]] <= q95) %>%
    ungroup() %>%
    ggplot(aes(x = sample, y = .data[[y_var]])) +
    geom_violin(aes(fill = X3243), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = c("firebrick", "darkslategray3", "goldenrod2")) +
    geom_point(
      data = . %>% filter(symbol == "Total" & barcode == pseudobulk_filter),
      aes(group = interaction(sample, X3243)), shape = 23, size = 1, fill = "white",
      position = position_dodge(width = 0.9)
    ) +
    theme_bw() +
    labs(x = NULL, y = y_label) +
    theme(
      aspect.ratio = 1, axis.title.x = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
  return(plot)
}

# Prepare data
filtered_data <- rbind(
  scmtMMB %>% filter(!sample %in% c("H05", "H47")),
  scmtMMB_ex,
  scmtMMB_exMax
)

# Generate P1 and save plots
P1 <- create_violin_plot(filtered_data, "mutation_per_MB", "Total scMPM", "pseudobulk", "../plot/22_scMPM_total_pm_MELAS_q5_q95.pdf")
P1.nl <- P1 + Seurat::NoLegend()
P1.leg <- ggpubr::get_legend(P1) %>% ggpubr::as_ggplot()
ggsave(plot = P1.nl, "../plot/22_scMPM_total_pm_MELAS_q5_q95.pdf", width = 3, height = 3)
ggsave(plot = P1.leg, "../plot/22_scMPM_total_pm_MELAS_leg_q5_q95.pdf", width = 3, height = 3)

# Generate P2 and save plots
P2 <- create_violin_plot(filtered_data, "MSS_weighted", "Total scwMSS", "pseudobulk", "../plot/22_scwMSS_total_pm_MELAS_q5_q95.pdf")
P2.nl <- P2 + Seurat::NoLegend()
P2.leg <- ggpubr::get_legend(P2) %>% ggpubr::as_ggplot()
ggsave(plot = P2.nl, "../plot/22_scwMSS_total_pm_MELAS_q5_q95.pdf", width = 3, height = 3)


theme_1_nogrid <- theme(aspect.ratio=1/1, axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

scmtMMB_exMax %>% dplyr::filter(symbol == "Total") %>%
  ggplot(aes(x = sample, y = MSS_weighted)) +
  #geom_jitter(aes(fill = sample), size = 0.1, color = "grey") + 
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  theme_bw() + xlab(NULL) + ylab("Total scwMSS")+
  theme_1_nogrid + Seurat::NoLegend()


scmtMMB_exMax %>% dplyr::filter(symbol == "Total") %>%
  ggplot(aes(x = sample, y = mutation_per_MB)) +
  geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
  scale_fill_manual(values = color_vec) +
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), 
             aes(group = sample), shape = 23, size = 1, fill = "white") + # overlay bulk data as point
  theme_bw() + xlab(NULL) + ylab("Total scMPM") +
  theme_1_nogrid+ Seurat::NoLegend() 
