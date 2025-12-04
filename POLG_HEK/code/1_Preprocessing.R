setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
library(Matrix)
library(Seurat)
library(Signac)
library(data.table)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(SummarizedExperiment)
library(ggplot2)
set.seed(123)
'%ni%' <- Negate('%in%')
suppressWarnings(annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

## Setup the Seurat Object
counts <-  Read10X_h5(filename = "../POLG_HEK_large_data_files/input/cellranger/filtered_peak_bc_matrix.h5")
fragpath <- "../POLG_HEK_large_data_files/input/cellranger/fragments.tsv.gz"
metadata <- read.csv(file="../POLG_HEK_large_data_files/input/cellranger/singlecell.csv", header = TRUE, row.names = 1)

# Create ChromatinAssay for ATAC data
atac_assay <- CreateChromatinAssay( 
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = fragpath,
  annotation = annotations)

# Create Seurat sbject
HEK <- CreateSeuratObject(counts = atac_assay,assay = 'ATAC', meta.data = metadata)
rm(counts, atac_assay, metadata)


## Quality Control
# nucleus genome ATAC QC
HEK$pct_reads_in_peaks <- HEK$peak_region_fragments / HEK$passed_filters * 100
HEK <- TSSEnrichment(HEK, fast = F)
HEK <- NucleosomeSignal(HEK)
HEK$log10_nCount_ATAC <- log10(HEK$nCount_ATAC)

# load mitochondrial data
mito.SE <- readRDS("../POLG_HEK_large_data_files/input/mgatk/mgatk.rds")[,Cells(HEK)]
HEK$mtDNA_depth <- mito.SE$depth
HEK$log10_mtDNA_depth <- log10(HEK$mtDNA_depth)
rm(mito.SE)

# visulize QC metrics 
VlnPlot(object = HEK,
        feature = c("TSS.enrichment",
                    "nucleosome_signal",
                    "pct_reads_in_peaks",
                    "nCount_ATAC",
                    "log10_mtDNA_depth"), pt.size = 0, ncol = 5)

## Demultiplexing with hashtag oligos (HTOs)
# Import HTO
import_hto<- function(){
  mtx <- fread(paste0("../POLG_HEK_large_data_files/input/featurecounts/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- Matrix::sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0("../POLG_HEK_large_data_files/input/featurecounts/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("../POLG_HEK_large_data_files/input/featurecounts/featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}


hto <- import_hto()
colnames(hto) <- paste0(colnames(hto), "-1")
HEK <- HEK[,Cells(HEK)%in% colnames(hto)]
hto <-hto[,Cells(HEK)]
HEK[["HTO"]] <- CreateAssayObject(hto[,Cells(HEK)])

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
HEK <- NormalizeData(HEK, assay = "HTO", normalization.method = "CLR")
HEK <- HTODemux(HEK, assay = "HTO", positive.quantile = 0.8)

# Global classification results
table(HEK$hash.ID)


# Visualize enrichment for selected HTOs with ridge plots
color_vec <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")
HEK$hash.ID <- factor(HEK$hash.ID, levels = c("CTRL", "KI36", "KIA2", "Negative","Doublet"))
Idents(HEK) <- "hash.ID"
p1 <- RidgePlot(HEK, assay = "HTO", features = rownames(HEK[["HTO"]])[1:3], ncol = 3, cols = c(color_vec, "Doublet" = "black", "Negative" = "grey"))
p1
ggsave(plot = p1, "../plot/Ext_Fig1b_11_HTO_demultiplex.pdf", width = 12, height = 4)

p <- HTOHeatmap(HEK, assay = "HTO") + scale_fill_viridis_c(limits = c(0,2))
     # scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "pink", 
     #                       high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 1, 
     #                       guide = "colourbar", aesthetics = "fill", limits = c(0,2))
p  # Ext_Fig1c

HEK <- readRDS("../POLG_HEK_large_data_files/output/1_POLG_HEK_preprocesssed_seurat.rds")

# visulize QC metrics 
VlnPlot(object = HEK,
        feature = c("TSS.enrichment",
                    "nucleosome_signal",
                    "pct_reads_in_peaks",
                    "log10_nCount_ATAC",
                    "log10_mtDNA_depth"),
        pt.size = 0, ncol = 5, cols = c(color_vec, "Doublet" = "black", "Negative" = "grey"))


## QC filtering
HEK_ori <- HEK
HEK <- subset(
  x = HEK,
  subset = nCount_ATAC > 1000 & nCount_ATAC < 150000 &
    pct_reads_in_peaks > 15 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2.5 &
    mtDNA_depth >= 5)

HEK

VlnPlot(object = HEK,
        feature = c("TSS.enrichment",
                    "nucleosome_signal",
                    "pct_reads_in_peaks",
                    "log10_nCount_ATAC",
                    "log10_mtDNA_depth"),
        pt.size = 0.01, ncol = 5, cols = c(color_vec, "Doublet" = "black", "Negative" = "grey"))

## Latent semantic indexing (LSI): Normalization (TF-IDF) followed by linear dimensional reduction (SVD)
set.seed(123)
DefaultAssay(HEK) <- "ATAC"
HEK <- RunTFIDF(HEK) # term frequency-inverse document frequency
HEK <- FindTopFeatures(HEK, min.cutoff = 'q0')
HEK <- RunSVD(HEK) # singular value decomposition

## Non-linear dimension reduction and clustering
set.seed(123)
HEK <- RunUMAP(object = HEK, reduction = 'lsi', dims = 2:30, reduction.name = "atac.umap", reduction.key = "atacUMAP_")
HEK <- FindNeighbors(object = HEK, reduction = 'lsi', dims = 2:30)
HEK <- FindClusters(object = HEK, verbose = FALSE, algorithm = 3)

DimPlot(object = HEK, group.by = "hash.ID", cols = c(color_vec, "Doublet" = "black", "Negative" = "grey")) 
DimPlot(object = HEK, label = TRUE) 

## Refine cell assignment
HEK <- subset(HEK, hash.ID != "Doublet" & hash.ID != "Negative")

HEK$condition <- data.frame(hash.ID = HEK$hash.ID, seurat_clusters = HEK$seurat_clusters) %>% 
  reframe(condition = case_when(HEK$hash.ID == "CTRL" & HEK$seurat_clusters %in% c(0,7,9) ~ "CTRL",
                                HEK$hash.ID == "KI36" & HEK$seurat_clusters %in% c(3,4,5,8) ~ "KI36",
                                HEK$hash.ID == "KIA2" & HEK$seurat_clusters %in% c(1,2) ~ "KIA2",
                                TRUE ~ "conflicted")) %>% pull(condition)

HEK <- subset(HEK, condition != "conflicted")
table(HEK$condition)
Idents(HEK) <- "condition"
levels(HEK) <- c("CTRL", "KI36", "KIA2")
table(Idents(HEK))
DimPlot(object = HEK, cols = color_vec) 

p2 <- DimPlot(object = HEK, cols = color_vec) + theme_void() + NoLegend()
ggsave(plot = p2, "../plot/Fig1b_12_UMAP.pdf", width = 5, height = 5)

VlnPlot(object = HEK,
        feature = c("TSS.enrichment",
                    "nucleosome_signal",
                    "pct_reads_in_peaks",
                    "nCount_ATAC",
                    "log10_mtDNA_depth"),
        pt.size = 0, ncol = 5, cols = color_vec)

p4 <- VlnPlot(object = HEK,
              feature = c("TSS.enrichment",
                          "nucleosome_signal",
                          "pct_reads_in_peaks",
                          "nCount_ATAC",
                          "log10_mtDNA_depth"),
              pt.size = 0, ncol = 5, cols = color_vec)

ggsave(plot = p4, "../plot/Ext_Fig1d_13_QC_filtered.pdf", width = 15, height = 4)


FeaturePlot(HEK, "CTRL") + scale_colour_viridis_c()
FeaturePlot(HEK, "KI36") + scale_colour_viridis_c()
FeaturePlot(HEK, "KIA2") + scale_colour_viridis_c()

# mtDNA_depth on UMAP

p4 <- FeaturePlot(object = HEK, features = "mtDNA_depth", min.cutoff = "q5", max.cutoff = "q95") + 
  scale_color_viridis_c() +
  theme_void() + ggtitle(NULL)
p4.nl <- p4 + theme_void() + Seurat::NoLegend() 
p4.leg <- ggpubr::get_legend(p4) %>% ggpubr::as_ggplot()
ggsave(plot = p4.nl, "../plot/Fig1c_14_mtDepth_UMAP_q5_q95.pdf", width = 5, height = 5)
ggsave(plot = p4.leg,"../plot/Fig1c_14_mtDepth_UMAP_leg.pdf", width = 3, height = 2)

# normalized mtDNA_depth on UMAP

HEK$mtDNA_depth_norm <- HEK$mtDNA_depth/HEK$nCount_ATAC
p4norm <- FeaturePlot(object = HEK, features = "mtDNA_depth_norm", min.cutoff = "q5", max.cutoff = "q95") + 
  scale_color_viridis_c() +  # min 5, max 200
  theme_void() + ggtitle(NULL)
p4norm.nl <- p4norm + theme_void() + Seurat::NoLegend() 
p4norm.leg <- ggpubr::get_legend(p4norm) %>% ggpubr::as_ggplot()
ggsave(plot = p4norm.nl, "../plot/Fig1d_14_norm_mtDepth_UMAP_q5_q95.pdf", width = 5, height = 5)
ggsave(plot = p4norm.leg,"../plot/Fig1d_14_norm_mtDepth_UMAP_leg.pdf", width = 3, height = 2)

# nCount_ATAC on UMAP
pA <- FeaturePlot(object = HEK, features = "nCount_ATAC", min.cutoff = "q5", max.cutoff = "q95") + 
  scale_color_viridis_c() +
  theme_void() + ggtitle(NULL)
pA.nl <- pA + theme_void() + Seurat::NoLegend() 
pA.leg <- ggpubr::get_legend(pA) %>% ggpubr::as_ggplot()
ggsave(plot = pA.nl, "../plot/Fig1c_14_nCount_ATAC_UMAP_q5_q95.pdf", width = 5, height = 5)
ggsave(plot = pA.leg,"../plot/Fig1c_14_nCount_ATAC_UMAP_leg.pdf", width = 3, height = 2)

p5 <- ggplot(data = HEK@meta.data, aes(x=nCount_ATAC, y = mtDNA_depth, color = condition)) +
  geom_point(size=0.01) + scale_color_manual(values = color_vec) + geom_smooth(method = "lm") +
  facet_wrap(~condition, scale = "free") +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0)) + 
  theme_classic() + NoLegend()

ggsave(plot = p5, "../plot/Ext_Fig1e_15_cor_ATAC_mtDepth.pdf", width = 6, height = 2)

p7 <- ggplot(data = HEK@meta.data, aes(x=condition, y = mtDNA_depth/nCount_ATAC)) +
  geom_violin(scale = "width", aes(fill = condition)) + geom_boxplot(width=0.1, outlier.size=0) + 
  scale_fill_manual(values = color_vec) + ylab("Normalized mtDNA depth") +
  theme_classic() + NoLegend() + scale_y_log10()

ggsave(plot = p7, "../plot/Fig1d_16_normalized_mtDNA_depth.pdf", width = 2.5, height = 2.5)

HEK@meta.data %>% group_by(condition) %>% summarise(median(mtDNA_depth/nCount_ATAC))
HEK@meta.data %>% group_by(condition) %>% summarise(sd(mtDNA_depth/nCount_ATAC))

# CTRL: 0.00258; KI36: 0.00756; KIA2: 0.00519 
# shallow: CTRL: 0.00265; KI36: 0.00732; KIA2: 0.00507 

library(rstatix)

stat.test <- HEK@meta.data %>% select(condition, mtDNA_depth_norm) %>% 
  wilcox_test(mtDNA_depth_norm ~ condition) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
stat.test
write.csv(stat.test, file = "../output/1_mtDNA_depth_norm_wilcox.test.csv", 
          quote = FALSE, row.names = TRUE)


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

feature.stat <- summarizer(HEK@meta.data, numeric_cols = c("nCount_ATAC", "mtDNA_depth","mtDNA_depth_norm"), condition) %>% 
  tidyr::pivot_longer(cols = -condition, names_to = "name", values_to = "value") %>% 
  mutate(
    feature = sub("^(.*)_[^_]+$", "\\1", name),
    stat = sub("^.*_([^_]+)$", "\\1", name) 
  ) %>% 
  select(condition, feature, stat, value) %>% 
  tidyr::pivot_wider(names_from = condition, values_from = value) %>% 
  mutate(All=unlist(summarizer(HEK@meta.data, numeric_cols = c("nCount_ATAC", "mtDNA_depth","mtDNA_depth_norm"))))
write.csv(feature.stat, file = "../output/1_features_stats.csv", 
          quote = FALSE, row.names = TRUE)

HEK <- AddMetaData(HEK, metadata = HEK@reductions$atac.umap@cell.embeddings)
saveRDS(object = HEK, file = "../POLG_HEK_large_data_files/output/1_POLG_HEK_preprocesssed_seurat.rds")
write.csv(x = HEK@meta.data, file = "../output/1_HEK_metadata.csv", row.names = TRUE, quote = FALSE)


# extract cell barcodes for each line
CTRL_CB <- WhichCells(HEK, idents = "CTRL")
KI36_CB <- WhichCells(HEK, idents = "KI36")
KIA2_CB <- WhichCells(HEK, idents = "KIA2")
CTRL <- subset(HEK, cells = CTRL_CB)
KI36 <- subset(HEK, cells = KI36_CB)
KIA2 <- subset(HEK, cells = KIA2_CB)

cor(CTRL$nCount_ATAC, y = CTRL$mtDNA_depth, method = "spearman") #0.4655705
cor(KI36$nCount_ATAC, y = KI36$mtDNA_depth, method = "spearman") #0.6427327
cor(KIA2$nCount_ATAC, y = KIA2$mtDNA_depth, method = "spearman") #0.6713671
