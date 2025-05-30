setwd("~/Ludwig_lab/scmtMMB/PBMC/code/")

## Load packages and genome annotation
library(Matrix)
# remotes::install_github("satijalab/seurat", "feat/dictionary", quiet = TRUE)
library(Seurat) 
library(Signac)
library(SeuratDisk)
library(data.table)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(SummarizedExperiment)
library(ggplot2)
set.seed(123)
'%ni%' <- Negate('%in%')
source("1_helper function.R")
source("../../global_func/variant_calling.R")
source("../../global_func/mtGene_colo.R")
color_vec <- c("M80" = "darkorchid4","M60" = "firebrick4", "M35" = "violetred3", "M29" = "firebrick1", "H47" = "blue3", "H05" = "dodgerblue3")

## Setup the Seurat Object and Quality control

IDs <- list.files("../PBMC_large_data_files/input/")

lapply(IDs, function(SID){
  pbmc <- load10XData(SID) %>% get_QC()
  mitoSE <- readRDS(paste0("../PBMC_large_data_files/input/", SID,"/mgatk/mgatk.rds"))[,Cells(pbmc)]
  pbmc$mtDNA_depth <- mitoSE$depth
})

suppressMessages({
  H05 <- load10XData("Healthy_Y05") %>% get_QC() ;H05$SID <- "H05"
  H47 <- load10XData("Healthy_Y47") %>% get_QC() ;H47$SID <- "H47"
  M29 <- load10XData("MELAS_Y29") %>% get_QC() ;M29$SID <- "M29"
  M35 <- load10XData("MELAS_Y35") %>% get_QC() ;M35$SID <- "M35"
  M60 <- load10XData("MELAS_Y60") %>% get_QC() ;M60$SID <- "M60"
  M80 <- load10XData("MELAS_Y80") %>% get_QC() ;M80$SID <- "M80"
  
  rm(obj.multi, obj.rna, obj.rna.ext)
})


H05.mitoSE <- readRDS("../PBMC_large_data_files/input/Healthy_Y05/mgatk/mgatk.rds")[,Cells(H05)];H05$mtDNA_depth <- H05.mitoSE$depth
H47.mitoSE <- readRDS("../PBMC_large_data_files/input/Healthy_Y47/mgatk/mgatk.rds")[,Cells(H47)];H47$mtDNA_depth <- H47.mitoSE$depth
M29.mitoSE <- readRDS("../PBMC_large_data_files/input/MELAS_Y29/mgatk/mgatk.rds")[,Cells(M29)];M29$mtDNA_depth <- M29.mitoSE$depth
M35.mitoSE <- readRDS("../PBMC_large_data_files/input/MELAS_Y35/mgatk/mgatk.rds")[,Cells(M35)];M35$mtDNA_depth <- M35.mitoSE$depth
M60.mitoSE <- readRDS("../PBMC_large_data_files/input/MELAS_Y60/mgatk/mgatk.rds")[,Cells(M60)];M60$mtDNA_depth <- M60.mitoSE$depth
M80.mitoSE <- readRDS("../PBMC_large_data_files/input/MELAS_Y80/mgatk/mgatk.rds")[,Cells(M80)];M80$mtDNA_depth <- M80.mitoSE$depth

# visulize QC metrics 
plotQC <- function(SO){
  p <- VlnPlot(object = SO, feature = c("TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks",
                                        "nCount_ATAC", "mtDNA_depth"), ncol = 5, pt.size = 0)
  # print(p)
  ggsave(plot = p, paste0("../plot/11_QC_",unique(SO$SID),".pdf"), height = 4, width = 15)
}

plotQC(H05); plotQC(H47); plotQC(M29); plotQC(M35); plotQC(M60); plotQC(M80)

QC_filter <- function(object){
  object <- subset(x = object, subset = nCount_ATAC > 1000 & nCount_ATAC < 20000 &
                     pct_reads_in_peaks > 25 & nucleosome_signal < 4 & TSS.enrichment > 2 & mtDNA_depth >= 10)
  return(object)
}

H05 <- QC_filter(H05); H47 <- QC_filter(H47); M29 <- QC_filter(M29); 
M35 <- QC_filter(M35); M60 <- QC_filter(M60); M80 <- QC_filter(M80)

plotQCf <- function(SO){
  p <- VlnPlot(object = SO, feature = c("TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks",
                                        "nCount_ATAC", "mtDNA_depth"), ncol = 5, pt.size = 0)
  ggsave(plot = p, paste0("../plot/12_QC_",unique(SO$SID),"filtered.pdf"), height = 4, width = 15)
}
plotQCf(H05); plotQCf(H47); plotQCf(M29); plotQCf(M35); plotQCf(M60); plotQCf(M80)

H05.mitoSE <- H05.mitoSE[,Cells(H05)]; H47.mitoSE <- H47.mitoSE[,Cells(H47)]; M29.mitoSE <- M29.mitoSE[,Cells(M29)]; 
M35.mitoSE <- M35.mitoSE[,Cells(M35)]; M60.mitoSE <- M60.mitoSE[,Cells(M60)]; M80.mitoSE <- M80.mitoSE[,Cells(M80)]

pull_coverage <- function(SE, cells, resolution = 5){
  zoo::rollmean(rowMeans(assays(SE)[['coverage']]), resolution)
}

cov_df <- data.frame(
  pos = zoo::rollmean(1:16569, 5),
  H05 = pull_coverage(H05.mitoSE), H47 = pull_coverage(H47.mitoSE), M29 = pull_coverage(M29.mitoSE),
  M35 = pull_coverage(M35.mitoSE), M60 = pull_coverage(M60.mitoSE), M80 = pull_coverage(M80.mitoSE)
)

write.csv(cov_df, "../output/1_mean_cov.csv", quote = FALSE, row.names = TRUE)

p3 <- cov_df %>% reshape2::melt(id.vars = "pos") %>% # dplyr::filter(value > 5) %>%
  ggplot() +
  geom_line(aes(x = pos, y = value, color = variable)) +
  scale_y_continuous(limits = c(-20, 65), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,16596, by = 5000), expand = c(0,0)) + 
  scale_color_manual(values = color_vec) +
  labs(x = "Position on mtDNA chromosome", y = "Roll mean coverage") + geom_rect(data = mitochondrial_data, 
                                                                                 aes(xmin = Starting, xmax = Ending, ymin = -3, ymax = 0, fill = color), alpha = 1) +
  scale_fill_manual(values = mtGene_colors) +
  coord_polar() + 
  guides(fill=guide_legend(ncol=2))+
  theme_minimal()

p3.nl <- p3 + Seurat::NoLegend()
p3.leg <- ggpubr::get_legend(p3) %>% ggpubr::as_ggplot()
ggsave("../plot/13_mtATAC_mtDNA_cov.pdf", p3.nl, width = 5, height = 5)
ggsave("../plot/13_mtATAC_mtDNA_cov_leg.pdf", p3.leg, width = 5, height = 5)

### mtDNA variant calling
H05.mmat <- call_mutations_mgatk(H05.mitoSE)
H47.mmat <- call_mutations_mgatk(H47.mitoSE)
M29.mmat <- call_mutations_mgatk(M29.mitoSE)
M35.mmat <- call_mutations_mgatk(M35.mitoSE)
M60.mmat <- call_mutations_mgatk(M60.mitoSE) 
M80.mmat <- call_mutations_mgatk(M80.mitoSE)

### mtDNA somatic variant
# Haplotypes were identify by selecting variant with VMR higher than 0.01 and high strand concordance (\> 0.65), 
# error-prone variants were excluded (301A\>C, 302A\>C, 310T\>C, and 316G\>C). 
# Top 10 high mean heteroplasmy variants in each cell line were visualized.

mgatk_filter <- function(mmat, n_cells = 1, vmr.thres = 0.01){
  var <- as.data.frame(rowData(mmat))
  var <- var %>% subset(n_cells_conf_detected >= n_cells &
                          strand_correlation > 0.65 &
                          vmr >= vmr.thres &
                          variant %ni% c("301A>C", "302A>C", "310T>C", "316G>C")) %>% pull(variant)
  mmat <- mmat[var,]
  return(mmat)
}

H05.mmat <- mgatk_filter(H05.mmat); H47.mmat <- mgatk_filter(H47.mmat); M29.mmat <- mgatk_filter(M29.mmat); 
M35.mmat <- mgatk_filter(M35.mmat); M60.mmat <- mgatk_filter(M60.mmat); M80.mmat <- mgatk_filter(M80.mmat)

saveData <- function(SO, mitoSE, mmat){
  md <- SO@meta.data %>% mutate(barcode = rownames(.)) %>% tibble::remove_rownames()
  write.table(x = md, file = paste0("../output/01_metadata_",unique(SO$SID),".tsv"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  saveRDS(SO,  paste0("../PBMC_large_data_files/output/01_seurat_",unique(SO$SID),".rds"))
  pos_cov <- assays(mitoSE)[['coverage']]
  saveRDS(pos_cov,  paste0("../PBMC_large_data_files/output/01_mtDNA_pos_coverage_",unique(SO$SID),".rds"))
  saveRDS(mmat,  paste0("../PBMC_large_data_files/output/01_mmat_",unique(SO$SID),".rds"))
}

saveData(H05, H05.mitoSE, H05.mmat)
saveData(H47, H47.mitoSE, H47.mmat)
saveData(M29, M29.mitoSE, M29.mmat)
saveData(M35, M35.mitoSE, M35.mmat)
saveData(M60, M60.mitoSE, M60.mmat)
saveData(M80, M80.mitoSE, M80.mmat)


