setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code")
library(SummarizedExperiment)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(Signac)
library(Seurat)
source("../../global_func/mtGene_color.R")


lapply(c("CTRL","KI36", "KIA2"), function(condt){
  set.seed(1234)
  mmtx <- readRDS(paste0("../POLG_HEK_large_data_files/output/2_mmtx.",condt,".rds"))
  
  VAF <- CreateSeuratObject(
    counts = assay(mmtx),
    assay = 'VAF'
  )
  
  VAF <- FindClonotypes(VAF)
  # DoHeatmap(VAF, features = VariableFeatures(VAF), slot = "data", disp.max = 0.1) +
  #  scale_fill_viridis_c()
  
  saveRDS(VAF, paste0("../output/scVAF/VAF_SO_", condt, ".rds"))
  
  CellID <- VAF@meta.data %>% tibble::rownames_to_column(var = "barcode") %>% 
    select(barcode, seurat_clusters) %>% group_by(seurat_clusters) %>%
    slice_head() %>% pull(barcode)
    write.csv(VAF@meta.data, paste0("../output/scVAF/", condt, ".csv"), quote = FALSE)
  
  VAF.long <- as.data.frame(assay(mmtx)[,CellID]) %>% 
    rename_with(.fn = ~paste0("Cell_", seq_along(.))) %>%
    rownames_to_column("variant") %>% 
    separate(variant, into = c("Position", "Change"), sep = "(?<=\\d)(?=\\D)") %>% 
    mutate(Position = as.numeric(Position)) %>% 
    pivot_longer(cols = -c("Position", "Change"), names_to = "ID", values_to = "VAF") 
  
  write.csv(VAF.long, paste0("../output/scVAF/", condt, "_coc.csv")) # cell of choice
  
  lapply(unique(VAF.long$ID),function(cell){
    
    filtered_data <- VAF.long %>% filter(ID == cell & VAF != 0)
    
    # Count the number of rows in the filtered dataframe
    num_rows <- nrow(filtered_data)
    
    p <- filtered_data %>% 
      ggplot() +
      geom_segment(aes(x=Position, xend=Position, y=0, yend=VAF)) +
      geom_point(aes(x = Position, y = VAF), size=0.1, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=2) +
      scale_x_continuous(limits = c(0,16569), expand = c(0,0)) +
      scale_y_continuous(limits = c(-0.5,1), expand = c(0,0)) +
      geom_hline(yintercept = c(0,1)) +
      geom_hline(yintercept = seq(0.1,0.9,0.1), linetype = "dashed",size= 0.2, color = "grey30") +
      geom_rect(data = mitochondrial_data, 
                aes(xmin = Starting, xmax = Ending, ymin = -0.2, ymax = -0.1, fill = color), alpha = 1) +
      scale_fill_manual(values = mtGene_colors) +
      coord_polar() + 
      guides(fill=guide_legend(ncol=2))+
      theme_void() + NoLegend() + theme(aspect.ratio=1/1) + 
      ggtitle(paste0(condt, "_",cell, "; number of variant:", num_rows))
      ggsave(plot = p, filename = paste0("../plot/scVAF/", condt,"/", cell,".pdf"), width = 5, height = 5)
  })
  
  count.long <- as.data.frame(assay(mmtx)[,CellID] * assays(mmtx)$coverage[,CellID]) %>% 
    rbind(depth = mmtx[,CellID]$depth,.) %>% 
    rename_with(.fn = ~paste0("Cell_", seq_along(.))) %>%
    rownames_to_column("variant") %>% 
    separate(variant, into = c("Position", "Change"), sep = "(?<=\\d)(?=\\D)") %>% 
    mutate(Position = as.numeric(Position)) %>% 
    pivot_longer(cols = -c("Position", "Change"), names_to = "ID", values_to = "count") 
  
  depth = count.long[!complete.cases(count.long),3:4]
  
  lapply(unique(count.long$ID),function(cell){
    
    
    p <- count.long %>% filter(ID == cell & count != 0) %>% 
      ggplot() +
      geom_segment(aes(x=Position, xend=Position, y=0, yend=count)) +
      geom_point(aes(x = Position, y = count), size=0.1, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=2) +
      scale_x_continuous(limits = c(0,16569), expand = c(0,0)) +
      scale_y_continuous(limits = c(-10, 150), expand = c(0,0)) +
      geom_hline(yintercept = c(0, as.numeric(depth[depth$ID == cell,"count"]))) +
      geom_hline(yintercept = seq(0,as.numeric(depth[depth$ID == cell,"count"]),10), linetype = "dashed",size= 0.2, color = "grey30") +
      geom_rect(data = mitochondrial_data, 
                aes(xmin = Starting, xmax = Ending, ymin = -3, ymax = -1, fill = color), alpha = 1) +
      scale_fill_manual(values = mtGene_colors) +
      coord_polar() + 
      guides(fill=guide_legend(ncol=2))+
      theme_void() + NoLegend() + theme(aspect.ratio=1/1) + 
      ggtitle(paste0(condt, "_",cell, ": mtDepth = ", as.numeric(depth[depth$ID == cell,"count"])))
      ggsave(plot = p, filename = paste0("../plot/scMut/", condt,"/", cell,".pdf"), width = 10, height = 10)
  })
  
})
