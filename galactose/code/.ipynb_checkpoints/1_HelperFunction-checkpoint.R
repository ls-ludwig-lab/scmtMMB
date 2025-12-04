import_hto<- function(duplicate){
  mtx <- fread(paste0("../data/",duplicate,"/featurecounts/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- Matrix::sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])

  rownames(matx) <- paste0(duplicate,"_", fread(paste0("../data/",duplicate,"/featurecounts/featurecounts.barcodes.txt"), header = FALSE)[[1]], "-1")
  colnames(matx) <- fread(paste0("../data/", duplicate,"/featurecounts/featurecounts.genes.txt"), header = FALSE)[[1]]
  return(t(matx))
}

import_hto<- function(duplicate){
  mtx <- fread(paste0("../data/",duplicate,"/featurecounts/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- Matrix::sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])

  rownames(matx) <- paste0(duplicate,"_", fread(paste0("../data/",duplicate,"/featurecounts/featurecounts.barcodes.txt"), header = FALSE)[[1]], "-1")
  colnames(matx) <- fread(paste0("../data/", duplicate,"/featurecounts/featurecounts.genes.txt"), header = FALSE)[[1]]
  return(t(matx))
}

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


VarPlot <- function(mmat, SampleID){
  counts <- as.data.frame(rowData(mmat)) %>% subset(n_cells_conf_detected >= 1 & 
                                                      strand_correlation > 0.65 & vmr > 0.005) %>% nrow()
  plt <- VariantPlot(as.data.frame(rowData(mmat)), vmr.threshold = 0.005) + 
    labs(title = SampleID, subtitle = paste0("Number of somatic variants: ", counts)) + 
    scale_y_log10(limits = c(-5,1), n.breaks = 7) + scale_x_continuous(limits = c(-0.6,1), n.breaks = 4)
  ggsave(plot = plt, paste0("../plot/Ext_9d_1_mtVarPlot_",SampleID,".pdf"), width = 3, height = 3)
  return(plt)
}