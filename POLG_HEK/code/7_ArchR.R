setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
suppressMessages({
  library(ArchR)
  library(Seurat)
  library(Signac)
  library(parallel)
  library(harmony)
  library(hexbin)
  library(MetBrewer)
  library(dplyr)
  library(ggrepel)
  source("ArchR_helperfunction.R")
  addArchRGenome("hg38")
  addArchRThreads(threads = 4) 
  color_vec <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")
})
set.seed(1234)

HEK.metadata <- read.csv("../output/1_HEK_metadata.csv", row.names = 1)

ArrowFiles <- createArrowFiles(
  inputFiles = "../POLG_HEK_large_data_files/input/cellranger/fragments.tsv.gz",
  sampleNames = "HEK", 
  minTSS = 3,
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

HEK.AR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "../POLG_HEK_large_data_files/output/ArchR/",
  copyArrows = FALSE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

idxSample <- BiocGenerics::which(HEK.AR$cellNames %in% paste0("HEK#",rownames(HEK.metadata))) 
cellsSample <- HEK.AR$cellNames[idxSample]
HEK.AR <- HEK.AR[cellsSample, ]

rownames(HEK.metadata) <- paste0("HEK#",rownames(HEK.metadata))
HEK.metadata <- HEK.metadata[cellsSample,]
HEK.AR$condition <- HEK.metadata$condition
HEK.AR$atacUMAP_1 <- HEK.metadata$atacUMAP_1
HEK.AR$atacUMAP_2 <- HEK.metadata$atacUMAP_2
HEK.AR$seurat_clusters <- HEK.metadata$seurat_clusters

HEK.AR$subclusters <- case_when(HEK.AR$seurat_clusters %in% c(0,7,9) ~ "CTRL",
                                HEK.AR$seurat_clusters %in% c(3,4) ~ "KI36_C1",
                                HEK.AR$seurat_clusters == 5 ~ "KI36_C3",
                                HEK.AR$seurat_clusters == 8 ~ "KI36_C2",
                                HEK.AR$seurat_clusters %in% c(1,2) ~ "KIA2",
                                TRUE ~ "conflicted") 

markerTest36 <- getMarkerFeatures(
  ArchRProj = HEK.AR, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "condition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "KI36", # higher in up-regulated
  bgdGroups = "CTRL" 
)

marker36 <- getMarkers(markerTest36, cutOff = "FDR <= 1 & abs(Log2FC) >= 0")
VP_36 <- plotVP(as.data.frame(marker36$KI36)) + 
  ggtitle("Differential Gene Activity (KI36 vs CTRL)") +
  scale_y_continuous(limits = c(NA, 45))
VP_36
ggsave(plot = VP_36, "../plot/Ext_Fig2a_VolcanoPlot_36.pdf", width = 8, height = 8, dpi = 300)

write.csv(as.data.frame(marker36$KI36), "../output/DGA/DA_36.csv", quote = FALSE, row.names = TRUE)

markerTestA2 <- getMarkerFeatures(
  ArchRProj = HEK.AR, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "condition",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "KIA2", # higher in up-regulated
  bgdGroups = "CTRL" 
)

markerA2 <- getMarkers(markerTestA2, cutOff = "FDR <= 1 & abs(Log2FC) >= 0")
VP_A2 <- plotVP(as.data.frame(markerA2$KIA2)) + ggtitle("Differential Gene Activity (KIA2 vs CTRL)")
VP_A2
ggsave(plot = VP_A2, "../plot/Ext_Fig2a_VolcanoPlot_A2.pdf", width = 8, height = 8, dpi = 300)

write.csv(as.data.frame(markerA2$KIA2), "../output/DGA/DA_A2.csv", quote = FALSE, row.names = TRUE)


marker36.sig.up <- getMarkers(markerTest36, cutOff = "FDR <= 0.05 & Log2FC >= 0.6")
marker36.sig.dn <- getMarkers(markerTest36, cutOff = "FDR <= 0.05 & Log2FC <= -0.6")
markerA2.sig.up <- getMarkers(markerTestA2, cutOff = "FDR <= 0.05 & Log2FC >= 0.6")
markerA2.sig.dn <- getMarkers(markerTestA2, cutOff = "FDR <= 0.05 & Log2FC <= -0.6")

writeLines(marker36.sig.up$KI36$name, "../output/DGA/marker36.sig.up.txt")
writeLines(marker36.sig.dn$KI36$name, "../output/DGA/marker36.sig.dn.txt")
writeLines(markerA2.sig.up$KIA2$name, "../output/DGA/markerA2.sig.up.txt")
writeLines(markerA2.sig.dn$KIA2$name, "../output/DGA/markerA2.sig.dn.txt")

markerUP <- intersect(marker36.sig.up$KI36$name, markerA2.sig.up$KIA2$name)
markerdn <- intersect(marker36.sig.dn$KI36$name, markerA2.sig.dn$KIA2$name)
writeLines(markerUP, "../output/DGA/marker.up.intersect.txt")
writeLines(markerdn, "../output/DGA/marker.dn.intersect.txt")

p <-  plotBrowserTrack(
  ArchRProj = HEK.AR, 
  groupBy = "condition", 
  geneSymbol = "MGME1", 
  upstream = 50000,
  downstream = 50000,
  pal = color_vec
)
grid::grid.newpage()
grid::grid.draw(p$MGME1)

plotPDF(plotList = p, 
        name = "MGME1-Tracks-Marker-Genes.pdf", # Ext_fig2b
        ArchRProj = HEK.AR, 
        addDOC = FALSE, width = 4, height = 4)

saveArchRProject(ArchRProj = HEK.AR, outputDirectory = "../POLG_HEK_large_data_files/output/ArchR/")

