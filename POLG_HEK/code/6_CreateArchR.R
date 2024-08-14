setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/POLG_HEK_large_data_files/output/ArchR")
suppressMessages({
  library(ArchR)
  library(Seurat)
  library(Signac)
  library(parallel)
  library(dplyr)
})

set.seed(1234)
addArchRGenome("hg38")
addArchRThreads(threads = 4) 
color_vec <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")

HEK.metadata <- read.csv("../../../output/1_HEK_metadata.csv", row.names = 1)

ArrowFiles <- createArrowFiles(
  inputFiles = "../../input/cellranger/fragments.tsv.gz",
  sampleNames = "HEK", 
  minTSS = 1, minFrags = 100, addTileMat = TRUE, addGeneScoreMat = TRUE
)

HEK.AR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = ".",
  copyArrows = FALSE)

idxSample <- BiocGenerics::which(HEK.AR$cellNames %in% paste0( "HEK#",rownames(HEK.metadata)))
cellsSample <- HEK.AR$cellNames[idxSample]
HEK.AR <- HEK.AR[cellsSample, ]
rownames(HEK.metadata) <- paste0("HEK#", rownames(HEK.metadata))
HEK.metadata <- HEK.metadata[cellsSample, ]
HEK.AR$Clusters <- HEK.metadata$seurat_clusters
HEK.AR$Condition <- HEK.metadata$condition
HEK.AR$CLusters2 <- ifelse(HEK.AR$Clusters %in% c(3,4), "KI36_1",
                           ifelse(HEK.AR$Clusters %in% c(1,6), "CTRL",
                                  ifelse(HEK.AR$Clusters %in% c(0,2), "KIA2",
                                         ifelse(HEK.AR$Clusters == "5", "KI36_2", "KI36_3"))))


# Transfer embedding from seurat
df <- data.frame(row.names=HEK.AR$cellNames, "seurat#UMAP1" = HEK.metadata$atacUMAP_1, "seurat#UMAP2" =  HEK.metadata$atacUMAP_2, check.names = FALSE)
HEK.AR@embeddings$seuratUMAP <- SimpleList(df = df, params = list())


saveArchRProject(ArchRProj = HEK.AR, outputDirectory = ".")

