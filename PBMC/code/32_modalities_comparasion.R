setwd("~/Ludwig_lab/scmtMMB/PBMC/code/")
library(Matrix)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(ggplot2)
library(viridis)
library(gridExtra)
require(scales)
'%ni%' <- Negate('%in%')
source("../../global_func/quantifyMMB.R")
source("../../global_func/variant_calling.R")
source("../../global_func/mtGene_color.R")

# Load the mitochondrial genome data from multiome
#ATAC_meta <- read.csv("../output/SCT_2306528004_Y63M_mtMultiome_pred.celltype.l2.csv", row.names = "X")
#ATAC_mgatk <- readRDS("../PBMC_large_data_files/input/SCT_2306528004_Y63M_mtMultiome/mgatk.rds")
ATAC_meta <- read.csv("../PBMC_large_data_files/input/Healthy_Y05_mtMultiome/Seurat/Hped_31697_mtMultiome_pred.celltype.l2.csv", row.names = "X")
ATAC_mgatk <- readRDS("../PBMC_large_data_files/input/Healthy_Y05_mtMultiome/mgatk/mgatk.rds")

ATAC_mmat <- call_mutations_mgatk(ATAC_mgatk)
ATAC_var <- as.data.frame(rowData(ATAC_mmat)) %>% subset(#n_cells_conf_detected >= 1 &
                                                         strand_correlation >= 0.65 & vmr >= 0.01 &
                                                         variant %ni% c("301A>C", "302A>C", "310T>C", "316G>C")) %>% pull(variant)
ATAC_mmat.filtered <- ATAC_mmat[ATAC_var, rownames(ATAC_meta)]
ATAC_cov_per_pos_cell <- assays(ATAC_mmat.filtered)[['coverage']]
ATAC.logic <- ATAC_cov_per_pos_cell >= 5
ATAC.mut.count <- assay(ATAC_mmat.filtered) * ATAC.logic

#RNA_mgatk <- readRDS("../PBMC_large_data_files/input/SCT_2306528004_Y63M_mtMultiome/mgatk.rna.rds")[,colnames(ATAC.mut.count)]
RNA_mgatk <- readRDS("../PBMC_large_data_files/input/Healthy_Y05_mtMultiome/mgatk/rna.mgatk.rds")[,colnames(ATAC.mut.count)]

RNA_mmat <- call_mutations_mgatk(RNA_mgatk)
RNA_mmat.filtered <- RNA_mmat[ATAC_var,]
RNA_cov_per_pos_cell <- assays(RNA_mmat.filtered)[['coverage']]
RNA.logic <- RNA_cov_per_pos_cell >= 5
RNA.mut.count <- assay(RNA_mmat.filtered) * RNA.logic

# visualize coverage

pull_coverage <- function(SE, cells, resolution = 5){
  zoo::rollmean(rowMeans(assays(SE)[['coverage']]), resolution)
}

cov_df <- data.frame(
  pos = zoo::rollmean(1:16569, 5),
  RNA = pull_coverage(RNA_mgatk),
  ATAC = pull_coverage(ATAC_mgatk))


p1 <- cov_df %>% reshape2::melt(id.vars = "pos") %>% # dplyr::filter(value > 5) %>%
  ggplot() +
  geom_line(aes(x = pos, y = value, color = variable)) +
  scale_y_continuous(limits = c(-20, 50), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,16596, by = 10000), expand = c(0,0)) + 
  scale_color_manual(values = c("red",  "blue")) +
  geom_hline(yintercept = c(0, 50)) +
  geom_hline(yintercept = seq(10,40,10), linetype = "dashed",size= 0.2, color = "grey30") +
  labs(x = "Position on mtDNA chromosome", y = "Roll mean coverage / cell")  + 
  geom_rect(data = mitochondrial_data, 
            aes(xmin = Starting, xmax = Ending, ymin = -3, ymax = 0, fill = color), alpha = 1) +
  scale_fill_manual(values = mtGene_colors) +
  coord_polar() + 
  theme_void() +
  guides(fill=guide_legend(ncol=2))+
  theme_minimal()

p2 <- p1 + Seurat::NoLegend()
leg <- ggpubr::get_legend(p1)
p3 <- ggpubr::as_ggplot(leg)

ggsave("../plot/32_H5_mtDNA_mtRNA_cov.pdf", p2, width = 5, height = 5)
ggsave("../plot/32_H5_mtDNA_mtRNA_cov_leg.pdf", p3, width = 5, height = 5)


# number of mutations
vaf.df <- data.frame(ATAC = colSums(ATAC.mut.count > 0),
                     RNA = colSums(RNA.mut.count > 0)) %>% 
  mutate(weights = ATAC/(ATAC + RNA + 0.01),
         log = weights > 0.5) #75%

p4 <- vaf.df %>% ggplot(aes(x = RNA, y = ATAC)) + geom_jitter(width = 0.2, height = 0.2, color = "grey") +
  geom_abline(slope = 1) + geom_density_2d(bins = 20) +
  xlim(-1, 17) + ylim(-1, 17) + theme_classic()

p5 <- vaf.df %>% 
  ggplot(aes(x = 1, y = weights))+
  geom_violin() +theme_classic()

ggsave("../plot/32_H5_mtDNA_mtRNA_mutation_counts.pdf", p4, width = 5, height = 5)
p4
ggsave("../plot/32_H5_mtDNA_mtRNA_mutation_weighs.pdf", p5, width = 1.5, height = 2)
