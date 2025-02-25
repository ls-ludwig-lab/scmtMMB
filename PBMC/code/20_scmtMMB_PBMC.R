setwd("~/Ludwig_lab/scmtMMB/PBMC/code/")
### Load packages
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
'%ni%' <- Negate('%in%'
                 )
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

p <- md %>% 
  ggplot(aes(x= refUMAP_1, y= refUMAP_2, color= predicted.celltype.l2)) +
  geom_point(size =0.05) + theme_void()
p.nl <- p + Seurat::NoLegend() 
p.leg <- ggpubr::get_legend(p) %>% ggpubr::as_ggplot()
ggsave(plot = p.nl, "../plot/2_refUMAP.pdf", width = 5, height = 5)
ggsave(plot = p.leg, "../plot/2_refUMAP.leg.pdf", width = 5, height = 5)


mtVar.bulk.df <- do.call(rbind, lapply(names(color_vec), function(SID){
  mmat <- readRDS(paste0("../PBMC_large_data_files/output/01_mmat_", SID,".rds"))
  rowData(mmat) %>% as_tibble() %>% mutate(sample = SID)
}))

mtVar.bulk.df <- mtVar.bulk.df %>% group_by(sample) %>% arrange(desc(mean)) %>% mutate(rank = 1:n()) 

mtVar.bulk.df %>% ggplot(aes(x = rank, y = mean)) + geom_point() + theme_classic() + facet_wrap(~sample) +
  ggtitle("Mean VAF")

scmtMMB <- do.call(rbind,lapply(names(color_vec),function(SID){
  mmat <- readRDS(paste0("../PBMC_large_data_files/output/01_mmat_", SID,".rds"))
  cov <- readRDS(paste0("../PBMC_large_data_files/output/01_mtDNA_pos_coverage_", SID,".rds"))
  rbind(
    quantify_MMB(mmat, cov, SID),
    quantify_MMB_psedobulk(mmat, cov, SID)
  )
}))


scmtMMB <- left_join(scmtMMB, md %>% select(predicted.celltype.l1, predicted.celltype.l2, 
                                            refUMAP_1, refUMAP_2, SID, mtDNA_depth, barcode, X3243A_G), 
                     by = c("barcode", "sample" = "SID"))

write.csv(scmtMMB, "../PBMC_large_data_files/output/2_scmtMMB.csv", quote = FALSE, row.names = TRUE)
#scmtMMB <- read.csv("../PBMC_large_data_files/output/2_scmtMMB.csv", row.names = "X")

scmtMMB_ex <- do.call(rbind,lapply(names(color_vec),function(SID){
  mmat <- readRDS(paste0("../PBMC_large_data_files/output/01_mmat_", SID,".rds"))
  mmat <- mmat[rownames(mmat) != "3243A>G", ]
  cov <- readRDS(paste0("../PBMC_large_data_files/output/01_mtDNA_pos_coverage_", SID,".rds"))
  rbind(
    quantify_MMB(mmat, cov, SID),
    quantify_MMB_psedobulk(mmat, cov, SID)
  )
}))

scmtMMB_ex <- left_join(scmtMMB_ex, md %>% select(predicted.celltype.l1, predicted.celltype.l2, 
                                                  refUMAP_1, refUMAP_2, SID, mtDNA_depth, barcode, X3243A_G), 
                        by = c("barcode", "sample" = "SID"))

write.csv(scmtMMB_ex, "../PBMC_large_data_files/output/2_scmtMMB_exclude_3243.csv", quote = FALSE, row.names = TRUE)
# scmtMMB_ex <- read.csv("../PBMC_large_data_files/output/2_scmtMMB_exclude_3243.csv", row.names = "X")
