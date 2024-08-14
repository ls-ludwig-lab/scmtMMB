setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
suppressMessages({
  library(ArchR)
  library(Seurat)
  library(Signac)
  library(parallel)
  library(dplyr)
  library(tidyr)
})
HEK.AR <- loadArchRProject(path = "../POLG_HEK_large_data_files/output/ArchR/")
#GSM <- getMatrixFromProject(HEK.AR)

markersGS <- getMarkerFeatures(
  ArchRProj = HEK.AR, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "CLusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerDF <- do.call(cbind,assays(markersGS)@listData) %>% 
  mutate(gene = rowData(markersGS)$name) %>% 
  pivot_longer(cols = 1:35, names_sep = "\\.", ,names_to = c( "variable", "condition"), values_to = "value") %>% 
  pivot_wider(id_cols = c("gene", "condition"), names_from = "variable", values_from = "value")

markerDF %>% filter(condition == "CTRL" & complete.cases(.)) %>% 
  mutate(color = ifelse(FDR < 0.05 & abs(Log2FC) >= 0.6, ifelse(Log2FC>= 0.6 ,'up','down'), 'no')) %>% 
  ggplot(aes(x= Log2FC, y = -log10(FDR), color = color)) + geom_point(size = 0.5) +
  scale_color_manual(values=c("slateblue", "grey","salmon")) + theme_bw() + NoLegend()

DGA_CvP <- markerDF %>% filter(condition == "CTRL" & complete.cases(.) & FDR < 0.05 & abs(Log2FC) >= 0.6) 
write.csv(DGA_CvP, file = "../output/7_DGA_CvP.csv",quote = FALSE, row.names = FALSE)

plotEmbedding(
  ArchRProj = HEK.AR, 
  colorBy = "GeneScoreMatrix", 
  name = "ZNF225", 
  embedding = "seuratUMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  #labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")


getDGA <- function(by = "Clusters", group1, group2){
  DGA <-getMarkerFeatures(
    ArchRProj = HEK.AR, 
    useMatrix = "GeneScoreMatrix",
    groupBy = by,
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2)
  
  df <- do.call(cbind, assays(DGA)@listData)
  colnames(df) <- names(assays(DGA))
  as.data.frame(cbind(rowData(DGA), df))
}

DGA_Cv363 <- getDGA(by = "CLusters2", "CTRL", "KI36_3")

DGA_Cv363 %>% mutate(color = ifelse(FDR < 0.05 & abs(Log2FC) >= 0.6, 
                          ifelse(Log2FC>= 0.6 ,'up','down'),
                          'no')) %>% 
  ggplot(aes(x= Log2FC, y = -log10(FDR), color = color)) + geom_point(size = 0.5) +
  scale_color_manual(values=c("slateblue", "grey","salmon")) + theme_bw() + NoLegend()

DGA_Cv362 <- getDGA(by = "CLusters2", "CTRL", "KI36_2")
DGA_Cv362 %>% mutate(color = ifelse(FDR < 0.05 & abs(Log2FC) >= 0.6, 
                                    ifelse(Log2FC>= 0.6 ,'up','down'),
                                    'no')) %>% 
  ggplot(aes(x= Log2FC, y = -log10(FDR), color = color)) + geom_point(size = 0.5) +
  scale_color_manual(values=c("slateblue", "grey","salmon")) + theme_bw() + NoLegend()

DGA_Cv361 <- getDGA(by = "CLusters2", "CTRL", "KI36_1")
DGA_Cv361 %>% mutate(color = ifelse(FDR < 0.05 & abs(Log2FC) >= 0.6, 
                                    ifelse(Log2FC>= 0.6 ,'up','down'),
                                    'no')) %>% 
  ggplot(aes(x= Log2FC, y = -log10(FDR), color = color)) + geom_point(size = 0.5) +
  scale_color_manual(values=c("slateblue", "grey","salmon")) + theme_bw() + NoLegend()

DGA_5v7 <- getDGA("5", "7")

DGA_5v7si <- DGA_5v7 %>% filter(FDR <= 0.05 & abs(Log2FC) >= 0.6)

write.csv(DGA_5v7si, file = "../output/7_DGA_5v7.csv",quote = FALSE, row.names = FALSE)
write.csv(DGA_5v7si, file = "../output/7_DGA_5v7si.csv",quote = FALSE, row.names = FALSE)

p1 <- DGA_5v7 %>% mutate(color = ifelse(FDR < 0.05 & abs(Log2FC) >= 0.6, 
                                  ifelse(Log2FC>= 0.6 ,'up','down'),
                                  'no')) %>% 
  ggplot(aes(x= Log2FC, y = -log10(FDR), color = color)) + geom_point(size = 0.5) +
  scale_color_manual(values=c("slateblue", "grey","salmon")) + theme_bw() + NoLegend()

ggsave(plot = p1, "../plot/71_DGA_cluster_5vs7.pdf", width = 3, height = 2)


markerDF %>% filter(condition == "KI36_3" & complete.cases(.)) %>% 
 ggplot(aes(x= Log2FC, y = -log10(FDR))) + geom_point()

