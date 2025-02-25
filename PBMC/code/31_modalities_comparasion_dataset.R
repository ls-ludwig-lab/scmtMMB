setwd("~/Ludwig_lab/scmtMMB/PBMC/code")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(ggplot2)
library(Matrix)
library(dplyr)
library(cowplot)
# other called packages SeuratDisk, details and DT
set.seed(123)
'%ni%' <- Negate('%in%')
suppressWarnings(annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"


# the 10x hdf5 file contains both data types. 
counts <- Read10X_h5("../PBMC_large_data_files/input/SCT_2306528004_Y63M_mtMultiome/filtered_feature_bc_matrix.h5")
fragpath <- "../PBMC_large_data_files/input/SCT_2306528004_Y63M_mtMultiome/atac_fragments.tsv.gz"

# load meta.data
metadata <- read.csv("../PBMC_large_data_files/input/SCT_2306528004_Y63M_mtMultiome/per_barcode_metrics.csv",
                     header = TRUE,
                     row.names = 1) 

# Create Seurat object
pbmc <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA", meta.data = metadata)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations)

pbmc


## Cellranger cell calling 
metadata <- metadata %>% 
  arrange(desc(atac_peak_region_fragments)) %>%
  mutate(rank_ATAC = 1:dim(metadata)[1],
         is_cell = ifelse(is_cell == 1 & excluded_reason == 0, "Cells", "Non-cells")) %>% 
  arrange(desc(gex_umis_count)) %>%
  mutate(rank_GEX = 1:dim(metadata)[1])

metadata %>%
  ggplot(aes(x = log10(rank_GEX), y = log10(gex_umis_count), col = factor(is_cell, levels = c("Non-cells", "Cells")))) + geom_line() +
  scale_y_continuous(breaks = 1:4, labels = c(10, 100, 1000, "10K"), expand = c(0, 0), limits = c(0, 5)) + 
  scale_x_continuous(breaks = c(2, 4, log10(20000), log10(50000), 5), labels = c(100, "10k", "20K","50K", "100K"), expand = c(0, 0), limits = c(0, 5.5)) + 
  scale_color_manual(values = c("blue2", "goldenrod2")) +
  theme_bw() + xlab("Barcodes") + ylab("GEX UMI count") +
  theme(legend.title=element_blank())


metadata %>%
  ggplot(aes(x = log10(rank_ATAC), y = log10(atac_peak_region_fragments), col = factor(is_cell, levels = c("Non-cells", "Cells")))) + geom_line() +
  scale_y_continuous(breaks = 1:4, labels = c(10, 100, 1000, "10K"), expand = c(0, 0), limits = c(0, 5)) + 
  scale_x_continuous(breaks = c(2, 4, log10(20000), log10(50000), 5), labels = c(100, "10k", "20","50", "100K"), expand = c(0, 0), limits = c(0, 5.5)) + 
  scale_color_manual(values = c("blue2", "goldenrod2")) +
  theme_bw() + xlab("Barcodes") + ylab("Fragments Overlapping Peaks") +
  theme(legend.title=element_blank())

metadata %>%
  ggplot(aes(x = atac_peak_region_fragments)) + 
  # geom_density(aes(y = stat(density) * 12000, col = factor(is_cell, levels = c("Non-cells", "Cells"))), alpha = 0.6, adjust = 3) +   
  geom_histogram(aes(fill = factor(is_cell, levels = c("Non-cells", "Cells"))),binwidth =0.05, alpha = 0.6, position = 'identity') +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = c("1", "10", 100, "1k", "10K", "100k"),
                limits = c(1, 1e5), expand = c(0,0)) + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = c("1", "10", 100, "1k", "10K", "100k"),
                limits = c(1, 1e5), expand = c(0,0)) +
  scale_fill_manual(values = c("blue2", "goldenrod2")) +
  scale_color_manual(values = c("blue2", "goldenrod2")) +
  theme_bw() + xlab("Fragments Per Barcode") + ylab("Barcodes") +
  theme(legend.title=element_blank())

## Quality Control
# nucleus genome ATAC QC
DefaultAssay(pbmc) <- "ATAC"
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

mtDepth <- read.table("../PBMC_large_data_files/input/SCT_2306528004_Y63M_mtMultiome/mgatk.depthTable.txt",
                      header = FALSE,
                      row.names = 1) 
pbmc$mtDNA_depth <- mtDepth[Cells(pbmc),]

# visulize QC metrics 
VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "mtDNA_depth", "percent.mt"),
  ncol = 3,
  pt.size = 0)

### QC filtering
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 25000 &
    nCount_RNA < 20000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 3 &
    percent.mt < 10
)

pbmc

#QC visualization after filtering**
VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "mtDNA_depth", "percent.mt"),
  ncol = 3,
  pt.size = 0)


## Gene expression and ATAC data processing
# GEX: SCTransform and PCA
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc) %>% RunUMAP(dims = 1:50, reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

# ATAC: latent semantic indexing (LSI)
DefaultAssay(pbmc) <- "ATAC"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc) %>% RunSVD() %>% 
  RunUMAP(reduction = 'lsi', dims = 2:50, reduction.name = "atac.umap", reduction.key = "atacUMAP_")


## Seurat v4 Reference Mapping
library(SeuratDisk)

# load PBMC reference

LabelTransfer <- function(SO){
  
  pbmc.cite <- LoadH5Seurat("../../../ref_data/pbmc_multimodal.h5seurat")
  
  DefaultAssay(SO) <- "SCT"
  
  # transfer cell type labels from reference to query
  anchors <- FindTransferAnchors(
    reference = pbmc.cite,
    query = SO,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  
  anchors <- FindTransferAnchors(
    reference = pbmc.cite,
    query = SO,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  
  mapped_SO <- MapQuery(
    anchorset = anchors,
    query = SO,
    reference = pbmc.cite,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2"),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
  
  SO$predicted.celltype.l1 <- mapped_SO$predicted.celltype.l1
  SO$predicted.celltype.l2 <- mapped_SO$predicted.celltype.l2
  # SO$predicted.celltype.l3 <- mapped_SO$predicted.celltype.l3
  SO <- AddMetaData(SO, mapped_SO@reductions$ref.umap@cell.embeddings, col.name = c("refUMAP_1", "refUMAP_2"))
  return(SO)
}

pbmc <- LabelTransfer(pbmc)

p <- FeatureScatter(pbmc,feature1 = "refUMAP_1", feature2 = "refUMAP_2", group.by = "predicted.celltype.l2", pt.size = 0.1) + 
  NoLegend() + ggtitle(NULL)
LabelClusters(p, id = "colors")


## Weighted Nearest Neighbor Analysis
pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)


#Visualize clustering based on each modality**
DimPlot(pbmc, reduction = "rna.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 2.5, repel = TRUE) + NoLegend() + ggtitle("RNA")
DimPlot(pbmc, reduction = "atac.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 2.5, repel = TRUE) + NoLegend() + ggtitle("ATAC")
DimPlot(pbmc, reduction = "wnn.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN") + NoLegend()

# Doublet removal 
## amulet
multiplets <- read.table("../PBMC_large_data_files/input/SCT_2306528004_Y63M_mtMultiome/amulet.txt") %>% t() %>% as.vector() #/MultipletCellIds_01.txt
multiplet.metadata <- rep("singlet", length(Cells(pbmc)))
names(multiplet.metadata) <- Cells(pbmc)
multiplet.metadata[multiplets] <- "multiplet"
pbmc$multiplet_amulet <- multiplet.metadata

DimPlot(pbmc, group.by = "multiplet_amulet", pt.size = 1, 
        cols = c("singlet" = "grey", "multiplet" = "purple"), 
        order = c("multiplet", "singlet"), reduction = "wnn.umap")

pbmc <- subset(pbmc, multiplet_amulet == "singlet") 
pbmc

saveRDS(pbmc, "../PBMC_large_data_files/output/SCT_2306528004_Y63M_mtMultiome_Seurat.rds")

celltype <- data.frame(predicted.celltype.l2 = pbmc$predicted.celltype.l2,
                         refUMAP_1 = pbmc$refUMAP_1,  refUMAP_2 = pbmc$refUMAP_2)
  
write.csv(celltype, "../output/SCT_2306528004_Y63M_mtMultiome_pred.celltype.l2.csv")

# Mitochondrial DNA variant calling
VlnPlot(pbmc, "mtDNA_depth", group.by = "predicted.celltype.l2", pt.size = 0.1) + NoLegend() + scale_y_log10() + geom_hline(yintercept=5)

#mtDNA coverage**
pbmc <- subset(pbmc, mtDNA_depth > 5)
mito <- readRDS("../PBMC_large_data_files/input/SCT_2306528004_Y63M_mtMultiome/mgatk.rds")
mito <- mito[,Cells(pbmc)]

pull_coverage <- function(SE, resolution = 5){
  zoo::rollmean(rowMeans(assays(SE)[['coverage']]), resolution)
}

cov_df <- data.frame(
  pos = zoo::rollmean(1:16569, 5),
  cov = pull_coverage(mito))

# visualize coverage
cov_df %>% reshape2::melt(id.vars = "pos") %>% # dplyr::filter(value > 5) %>%
  ggplot(aes(x = pos, y = value)) +
  geom_line() + theme_classic() +
  scale_x_continuous(breaks = seq(0,16596, by = 2500), expand = c(0,0)) +
  labs(x = "Position on mtDNA chromosome", y = "Roll mean coverage") 

source("../../global_func/variant_calling.R")
suppressWarnings(mmat <- call_mutations_mgatk(mito))
VariantPlot(
  as.data.frame(rowData(mmat)),
  min.cells = 1,
  concordance.threshold = 0.65,
  vmr.threshold = 0.01
)

var <- as.data.frame(rowData(mmat))
var <- var %>% subset(n_cells_conf_detected >= 1 &
                        strand_correlation > 0.65 &
                        # vmr >= vmr.thres &
                        variant %ni% c("301A>C", "302A>C", "310T>C", "316G>C")) %>% pull(variant)
mmat <- mmat[var,]

saveRDS(mmat, "../PBMC_large_data_files/output/SCT_2306528004_Y63M_mtMultiome_mmat.rds")



