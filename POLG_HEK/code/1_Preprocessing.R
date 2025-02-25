setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code")

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
# refallele <- data.frame(pos = 1:16569, as.data.frame(rowData(mito.SE)))
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
hto <-hto[,Cells(HEK)]
HEK[["HTO"]] <- CreateAssayObject(hto[,Cells(HEK)])

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
HEK <- NormalizeData(HEK, assay = "HTO", normalization.method = "CLR")
HEK <- HTODemux(HEK, assay = "HTO", positive.quantile = 0.8)

# Global classification results
table(HEK$hash.ID)
      # Doublet     KIA2     CTRL     KI36 Negative 
      # 1432     2084     1747     1810      346 

# Visualize enrichment for selected HTOs with ridge plots
color_vec <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")
HEK$hash.ID <- factor(HEK$hash.ID, levels = c("CTRL", "KI36", "KIA2", "Negative","Doublet"))
Idents(HEK) <- "hash.ID"
p1 <- RidgePlot(HEK, assay = "HTO", features = rownames(HEK[["HTO"]])[1:3], ncol = 3, cols = c(color_vec, "Doublet" = "black", "Negative" = "grey"))
p1
ggsave(plot = p1, "../plot/11_HTO_demultiplex.pdf", width = 12, height = 4)

# visulize QC metrics 
VlnPlot(object = HEK,
        feature = c("TSS.enrichment",
                    "nucleosome_signal",
                    "pct_reads_in_peaks",
                    "nCount_ATAC",
                    "log10_mtDNA_depth"),
        pt.size = 0, ncol = 5, cols = c(color_vec, "Doublet" = "black", "Negative" = "grey"))


## QC filtering
HEK <- subset(
  x = HEK,
  subset = nCount_ATAC > 1000 & nCount_ATAC < 20000 &
    pct_reads_in_peaks > 15 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2.5 &
    mtDNA_depth > 5)

HEK

VlnPlot(object = HEK,
        feature = c("TSS.enrichment",
                    "nucleosome_signal",
                    "pct_reads_in_peaks",
                    "nCount_ATAC",
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
  reframe(condition = case_when(HEK$hash.ID == "CTRL" & HEK$seurat_clusters %in% c(1,6) ~ "CTRL",
                                HEK$hash.ID == "KI36" & HEK$seurat_clusters %in% c(3,4,5,7) ~ "KI36",
                                HEK$hash.ID == "KIA2" & HEK$seurat_clusters %in% c(0,2) ~ "KIA2",
                                TRUE ~ "conflicted")) %>% pull(condition)

HEK <- subset(HEK, condition != "conflicted")
Idents(HEK) <- "condition"
levels(HEK) <- c("CTRL", "KI36", "KIA2")
table(Idents(HEK))
DimPlot(object = HEK, cols = color_vec) 

p2 <- DimPlot(object = HEK, cols = color_vec) + theme_void() + NoLegend()
ggsave(plot = p2, "../plot/12_UMAP.pdf", width = 5, height = 5)

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

ggsave(plot = p4, "../plot/13_QC_filtered.pdf", width = 15, height = 4)


# mtDNA_depth on UMAP

p4 <- FeaturePlot(object = HEK, features = "mtDNA_depth", min.cutoff = "q5", max.cutoff = "q95") + 
  scale_color_viridis_c() +
  theme_void() + ggtitle(NULL)
p4.nl <- p4 + theme_void() + Seurat::NoLegend() 
p4.leg <- ggpubr::get_legend(p4) %>% ggpubr::as_ggplot()
ggsave(plot = p4.nl, "../plot/14_mtDepth_UMAP_q5_q95.pdf", width = 5, height = 5)
ggsave(plot = p4.leg,"../plot/14_mtDepth_UMAP_leg.pdf", width = 3, height = 2)

# normalized mtDNA_depth on UMAP

HEK$mtDNA_depth_norm <- HEK$mtDNA_depth/HEK$nCount_ATAC
p4norm <- FeaturePlot(object = HEK, features = "mtDNA_depth_norm", min.cutoff = "q5", max.cutoff = "q95") + 
  scale_color_viridis_c() +  # min 5, max 200
  theme_void() + ggtitle(NULL)
p4norm.nl <- p4norm + theme_void() + Seurat::NoLegend() 
p4norm.leg <- ggpubr::get_legend(p4norm) %>% ggpubr::as_ggplot()
ggsave(plot = p4norm.nl, "../plot/14_norm_mtDepth_UMAP_q5_q95.pdf", width = 5, height = 5)
ggsave(plot = p4norm.leg,"../plot/14_norm_mtDepth_UMAP_leg.pdf", width = 3, height = 2)

# nCount_ATAC on UMAP
pA <- FeaturePlot(object = HEK, features = "nCount_ATAC", min.cutoff = "q5", max.cutoff = "q95") + 
  scale_color_viridis_c() +
  theme_void() + ggtitle(NULL)
pA.nl <- pA + theme_void() + Seurat::NoLegend() 
pA.leg <- ggpubr::get_legend(pA) %>% ggpubr::as_ggplot()
ggsave(plot = pA.nl, "../plot/14_nCount_ATAC_UMAP_q5_q95.pdf", width = 5, height = 5)
ggsave(plot = pA.leg,"../plot/14_nCount_ATAC_UMAP_leg.pdf", width = 3, height = 2)

p5 <- ggplot(data = HEK@meta.data, aes(x=nCount_ATAC, y = mtDNA_depth, color = condition)) +
  geom_point(size=0.01) + scale_color_manual(values = color_vec) + geom_smooth(method = "lm") +
  facet_wrap(~condition) +
  scale_x_continuous(limits = c(0, NA), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 150), expand = c(0,0)) + 
  theme_classic() + NoLegend()

ggsave(plot = p5, "../plot/15_cor_ATAC_mtDepth.pdf", width = 6, height = 2)

p7 <- ggplot(data = HEK@meta.data, aes(x=condition, y = mtDNA_depth/nCount_ATAC)) +
  geom_violin(scale = "width", aes(fill = condition)) + geom_boxplot(width=0.1, outlier.size=0) + 
  scale_fill_manual(values = color_vec) + ylab("Normalized mtDNA depth") +
  theme_classic() + NoLegend() + scale_y_log10()

ggsave(plot = p7, "../plot/16_normalized_mtDNA_depth.pdf", width = 2.5, height = 2.5)

HEK@meta.data %>% group_by(condition) %>% summarise(median(mtDNA_depth/nCount_ATAC))
# CTRL: 0.00298; KI36: 0.00741; KIA2: 0.00560

library(rstatix)

stat.test <- HEK@meta.data %>% select(condition, mtDNA_depth_norm) %>% 
  wilcox_test(mtDNA_depth_norm ~ condition) %>%
  adjust_pvalue() %>%
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



