setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
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
  scale_y_continuous(limits = c(-10, NA), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,16596, by = 5000), expand = c(0,0)) + 
  scale_color_manual(values = color_vec) +
  labs(x = "Position on mtDNA chromosome", y = "Roll mean coverage") + geom_rect(data = mitochondrial_data, 
                                                                                 aes(xmin = Starting, xmax = Ending, ymin = -10, ymax = 0, fill = color), alpha = 1) +
  scale_fill_manual(values = mtGene_colors) +
  # coord_polar() + 
  guides(fill=guide_legend(ncol=2))+
  theme_classic()

p1.nl <- p1 + Seurat::NoLegend()
leg <- ggpubr::get_legend(p1) %>% ggpubr::as_ggplot()
ggsave("../plot/Fig1e_21_mtDNA_cov.pdf", p1.nl, width = 5, height = 3)
ggsave("../plot/Fig1e_21_mtDNA_cov_leg.pdf", leg, width = 5, height = 5)

## Mitochondrial variant calling
mmat.CTRL <- call_mutations_mgatk(mitoSE.CTRL)
mmat.KI36 <- call_mutations_mgatk(mitoSE.KI36)
mmat.KIA2 <- call_mutations_mgatk(mitoSE.KIA2)

## Homoplasmic variant
mmat.CTRL 
mmat.KI36 
mmat.KIA2 

homoVAR <- function(mmat){
  as.data.frame(rowData(mmat)) %>% subset(n_cells_conf_detected >= 1 & 
                                                      strand_correlation > 0.65 & vmr <= 0.003) %>% pull(variant)
}

VarPlot <- function(mmat, SampleID){
  counts <- as.data.frame(rowData(mmat)) %>% subset(n_cells_conf_detected >= 1 & 
                                                      strand_correlation > 0.65 & vmr >= 0.003) %>% nrow()
  plt <- VariantPlot(as.data.frame(rowData(mmat)), vmr.threshold = 0.003) + 
    labs(title = SampleID, subtitle = paste0("Number of somatic variants: ", counts)) + 
    scale_y_log10(limits = c(-5,1), n.breaks = 7) + scale_x_continuous(limits = c(-0.6,1), n.breaks = 4)
  ggsave(plot = plt, paste0("../plot/Fig1f_22_VMR_",SampleID,".pdf"), width = 3, height = 3)
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
                         ifelse(vmr > 0.003 & mean < 0.9, "high-confident", "germline"), "error-prone"))
  table(var$type)
}

sumVar(mmat.CTRL)
# error-prone       germline high-confident 
# 1349             28            620 

sumVar(mmat.KI36)
# error-prone      germline high-confident 
# 11744             28           9655 

sumVar(mmat.KIA2)
#error-prone       germline high-confident 
#11181             27          11407 

mtMut <- function(mmat, n_cells = 1, vmr.thres = 0.003){
  var <- as.data.frame(rowData(mmat))
  var <- var %>% subset(n_cells_conf_detected >= n_cells &
                          strand_correlation > 0.65 & mean < 0.9 &
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




