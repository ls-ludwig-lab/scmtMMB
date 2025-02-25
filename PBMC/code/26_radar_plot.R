#09 calculate the enrichment of a variant within a specific cell type or lineage and to verify kurskal-wallis test

#phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)

#x, q vector of quantiles representing the number of white balls drawn
#without replacement from an urn which contains both black and white
#balls. # number of cells with mutation in the particular celltype 
#add -1 and lower.tail = FALSE because lower.tail: logical; if TRUE (default), probabilities are P[X <= x], 
#otherwise, P[X > x]

#m the number of white balls in the urn. #number of cells with the mutation

#n the number of black balls in the urn. #number of cells w/o the mutation

#k the number of balls drawn from the urn. # number of cells in the celltype
#https://stackoverflow.com/questions/8382806/hypergeometric-test-phyper

setwd("~/Ludwig_lab/scmtMMB/PBMC/code/")

#load libraries
library(dplyr)
library(stringr)
library(plyr)
library(ggplot2)
library(ggrepel)
library(fmsb)
library(SummarizedExperiment)
library(Matrix)

var_artifactX <- c('X301A.C', 'X302A.C', 'X309C.T','X310T.C', 'X316G.C', 'X3109T.C')
'%ni%' <- Negate('%in%')

scmtMMB <- read.csv("../PBMC_large_data_files/output/2_scmtMMB.csv", row.names = "X") %>% filter(barcode != "pseudobulk" & symbol == "Total") %>%
  mutate(
    celltype = case_when(
      predicted.celltype.l2 %in% c("CD8 TEM", "CD8 TCM") ~ "CD8 Eff/Mem",
      predicted.celltype.l2 %in% c("CD4 TEM", "CD4 TCM", "CD4 CTL") ~ "CD4 Eff/Mem",
      predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory", "Plasmablast") ~ "B cell",
      predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC") ~ "Monocyte/DC",
      grepl("NK|ILC", predicted.celltype.l2) ~ "NK/ILC",
      predicted.celltype.l2 %in% c("Doublet", "Eryth", "Platelet", "ASDC", NA) ~ "discard",
      TRUE ~ predicted.celltype.l2
    )
  ) %>% filter(celltype != "discard")

H5_mmat <- readRDS("../PBMC_large_data_files/output/01_mmat_H05.rds")
df_H05 <- scmtMMB %>% filter(sample == "H05") %>% mutate(VAF = assay(H5_mmat)["10599G>A",.$barcode]) 

M60_mmat <- readRDS("../PBMC_large_data_files/output/01_mmat_M60.rds")
df_M60 <- scmtMMB %>% filter(sample == "M60") %>% mutate(VAF = assay(M60_mmat)["10270T>C",.$barcode])

M80_mmat <- readRDS("../PBMC_large_data_files/output/01_mmat_M80.rds")
df_M80 <- scmtMMB %>% filter(sample == "M80") %>% mutate(VAF = assay(M80_mmat)["10398A>G",.$barcode])


df.list = list(df_H05, df_M60, df_M80)
df_meta = data.frame(
  name = c("H05", "M60", "M80"),
  var = c("10599G>A", "10270T>C", "10398A>G")
)



df.list[[3]] %>% 
  filter(symbol == "Total" & celltype %in% c("CD4 Naive", "CD4 Eff/Mem", "CD8 Naive", "CD8 Eff/Mem", "B cell", "Monocyte/DC", "NK/ILC")) %>% 
  ggplot(aes(x = celltype, y = VAF)) +
  geom_violin(scale = "width")

df.list[[2]] %>% 
  filter(symbol == "Total" & celltype %in% c("CD4 Naive", "CD4 Eff/Mem", "CD8 Naive", "CD8 Eff/Mem", "B cell", "Monocyte/DC", "NK/ILC")) %>% 
  ggplot(aes(VAF, color = celltype)) +
  stat_ecdf(geom = "step")

df.list[[3]] %>% 
  filter(symbol == "Total" & celltype %in% c("CD4 Naive", "CD4 Eff/Mem", "CD8 Naive", "CD8 Eff/Mem", "B cell", "Monocyte/DC", "NK/ILC")) %>% 
  ggplot(aes(VAF, color = celltype)) +
  stat_ecdf(geom = "step")

lapply(1:3, function(x){
  #load data
  celltype_mito_combined <- df.list[[x]] #combination of output of summarized experiment with celltype annotation, rownames = cellbc, colnames = mtDNA + celltype annotation
  
  ####################################
  #####----CELLTYPES LEVEL 1----######
  ####################################
  
  #make dataframe
  
  #enrichment analysis using phyper
  bias <- do.call(rbind, 
                  lapply(unique(celltype_mito_combined$celltype), function(ct){
                    q <- sum(celltype_mito_combined[celltype_mito_combined$celltype == ct,"VAF"]>0) 
                    m <- sum(celltype_mito_combined[,"VAF"]>0)
                    n <- nrow(celltype_mito_combined)-m
                    k <- nrow(celltype_mito_combined[celltype_mito_combined$celltype == ct,])
                    df <- data.frame(celltype = ct,bias = phyper(q, m, n, k, lower.tail = FALSE))
                  }
                  )    
  )
  
  #adj. p value!!
  bias$pvalue_adj <- p.adjust(bias$bias, method = "BH")
  
  write.csv(bias, paste0("../output/26_", df_meta[x,1],"_",df_meta[x,2], "_enrichment_analysis_phyper.csv"))
  
  # write.csv(bias, paste0(id, "enrichment_analysis_celltypes_level_1_phyper.csv"))
  
  #filter for the significant ones below 0.05 adj. p-value
  bias_filter <- bias[bias$pvalue_adj < 0.05,] %>% arrange(pvalue_adj) 
  
  #----------- make radial plots with all timepoints ----------------------------#
  
  #subset dataframe to do radial plots with these variants
  
  #check which cells have the mutation at all
  above_0 <- celltype_mito_combined[,c("celltype", "VAF")] %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarize(n = n(),
                     across(-n, ~sum( . > 0)))%>%
    ungroup()
  
  #check which cells have the mutation above 10 % heteroplasmy
  above_10 <- celltype_mito_combined[,c("celltype", "VAF")] %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarize(n = n(),
                     across(-n, ~sum( . > 0.1)))%>%
    ungroup()
  
  #check which cells have the mutation above 60 % heteroplasmy
  above_60 <- celltype_mito_combined[,c("celltype", "VAF")] %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarize(n = n(),
                     across(-n, ~sum( . > 0.6)))%>%
    ungroup()
  
  
  #take relative values
  above_0_2 <- above_0[,-1] #get rid off the predicted cell id for easier calculations
  above_10_2 <- above_10[,-1]
  above_60_2 <- above_60[,-1]
  
  #get relative values by dividing the whole dataframe by the number of observations
  relative_values_above_0 <- above_0_2/t(above_0_2[,"n"])
  relative_values_above_10 <- above_10_2/t(above_10_2[,"n"])
  relative_values_above_60 <- above_60_2/t(above_60_2[,"n"])
  
  #have the celltype annotation again
  rownames(relative_values_above_0) <- above_0$celltype
  rownames(relative_values_above_10) <- above_10$celltype
  rownames(relative_values_above_60) <- above_60$celltype
  
  
  
  
    #plot for each mutation the radials plot by making a graph for each
    test_radial_plot_above_0 <- as.data.frame(relative_values_above_0[,"VAF"])
    test_radial_plot_above_10 <- as.data.frame(relative_values_above_10[,"VAF"])
    test_radial_plot_above_60 <- as.data.frame(relative_values_above_60[,"VAF"])
    
    rownames(test_radial_plot_above_0) <- rownames(relative_values_above_0)
    rownames(test_radial_plot_above_10) <- rownames(relative_values_above_10)
    rownames(test_radial_plot_above_60) <- rownames(relative_values_above_60)
    #test_radial_plot$predicted.celltype.l1 <- NULL
    
    
    data <- as.data.frame(rbind(rep(max(t(test_radial_plot_above_0)),length(unique(celltype_mito_combined$celltype))),
                                rep(0,length(unique(celltype_mito_combined$celltype))), t(test_radial_plot_above_0)))
                                # t(test_radial_plot_above_0), t(test_radial_plot_above_10), t(test_radial_plot_above_60)))
    write.csv(data, paste0("../output/26_", df_meta[x,1],"_",df_meta[x,2], "_radial_prop.csv"))
    
    subdata <- data %>% select(`B cell`, `CD4 Eff/Mem`, `CD4 Naive`, `CD8 Eff/Mem`, `CD8 Naive`, `Monocyte/DC`, `NK/ILC`) 
    sub10 <- as.data.frame(t(test_radial_plot_above_0)) %>% 
              select(`B cell`, `CD4 Eff/Mem`, `CD4 Naive`, `CD8 Eff/Mem`, `CD8 Naive`, `Monocyte/DC`, `NK/ILC`)
    
    pdf(paste0("../plot/26_", df_meta[x,1],"_",df_meta[x,2], "_radial_plot.pdf"), width = 4, height = 4)
    
    radarchart(subdata, 
               cglty = 1,       
               cglcol = "gray", 
               cglwd = 1,       
               pcol= "red", pty =32, plwd = 4, 
               title = paste0(df_meta[x,1],"_",df_meta[x,2],", max_prop:",round(max(sub10),3)))
    dev.off()

   })
  
  