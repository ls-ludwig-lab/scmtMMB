library(data.table)
library(dplyr)
mtGene <- fread("~/Ludwig_lab/scmtMMB/database/Coding_region.tsv") %>% mutate(Length = Ending - Starting + 1)
mtGene$Positions <- sapply(1:37, function(i) mtGene$Starting[i]:mtGene$Ending[i])
mtRegion <- sapply(unique(mtGene$Region), function(x) unique(unlist(mtGene$Positions[mtGene$Region == x])))
MLC_score <- fread("~/Ludwig_lab/scmtMMB//database/MLC_score.tsv") # Mitochondrial local constraint (MLC) 


# quantify_MMB
quantify_MMB <- function(mmat, cov_per_pos_cell, SampleID){
  CellID <- colnames(mmat)
  var_meta <- left_join(as.data.frame(rowData(mmat)), MLC_score, by = c("position" = "Position"))
  
  # Calculate per position counts and MLC score
  mtCov_per_pos_cell <- cov_per_pos_cell[,CellID]
  mtMut_counts_per_pos_cell <- assays(mmat)$allele_frequency*assays(mmat)$coverage
  MLC_per_pos_cell <- sweep(assays(mmat)$allele_frequency > 0, 1, var_meta$MLC_pos_score, "*") 
  MLC_per_pos_cell_weighted <- sweep(assays(mmat)$allele_frequency, 1, var_meta$MLC_pos_score, "*") 
  
  # Calculate per cell counts, MMB and MLC score
  mtCov_cell <- colSums(mtCov_per_pos_cell)
  mtCopy_cell <- mtCov_cell/16569
  mtMut_counts_cell <- colSums(mtMut_counts_per_pos_cell)
  Mut_MB_cell <- mtMut_counts_cell/mtCov_cell * 1e6
  Mut_copy_cell <- mtMut_counts_cell/mtCopy_cell
  MSS <- colSums(MLC_per_pos_cell) # MLC score sum (MSS)
  MSS_weighted <- colSums(MLC_per_pos_cell_weighted) # MLC score sum (MSS)
  
  df1 <- data.frame(barcode = CellID, sample = SampleID, symbol = "Total", 
                    coverage = mtCov_cell, copy_number_abs = mtCopy_cell, mutation_counts = mtMut_counts_cell,
                    mutation_per_MB = Mut_MB_cell, mutation_per_copy = Mut_copy_cell, 
                    MSS = MSS, MSS_oe = MSS/sum(MLC_score$MLC_pos_score), MSS_weighted = MSS_weighted)
  
  calculate_gene_metrics <- function(i, region, symbol) {
    row_idx <- var_meta$position %in% region
    mtCov_gene_cell <- colSums(mtCov_per_pos_cell[region, ])
    mtCopy_gene_cell <- mtCov_gene_cell / length(region)
    
    if (sum(row_idx) > 1) {
      mtMut_counts_gene_cell <- colSums(mtMut_counts_per_pos_cell[row_idx, ])
      MSS_gene_cell <- colSums(MLC_per_pos_cell[row_idx, ])
      MSS_gene_cell_weighted <- colSums(MLC_per_pos_cell_weighted[row_idx, ])
    } else if (sum(row_idx) == 1) {
      mtMut_counts_gene_cell <- mtMut_counts_per_pos_cell[row_idx, ]
      MSS_gene_cell <- MLC_per_pos_cell[row_idx, ]
      MSS_gene_cell_weighted <- MLC_per_pos_cell_weighted[row_idx, ]
    } else {
      mtMut_counts_gene_cell <- MSS_gene_cell <- MSS_gene_cell_weighted <- rep(0, length(CellID))
    }
    
    if (sum(row_idx) > 0) {
      MSS_oe <- MSS_gene_cell / sum(MLC_score$MLC_pos_score[MLC_score$Position %in% region])
    } else {
      MSS_oe <- 0
    }
    
    Mut_MB_cell <- mtMut_counts_gene_cell / mtCov_gene_cell * 1e6
    Mut_copy_cell <- mtMut_counts_gene_cell / mtCopy_gene_cell
    
    data.frame(
      barcode = CellID, sample = SampleID, symbol = symbol, 
      coverage = mtCov_gene_cell, copy_number_abs = mtCopy_gene_cell, 
      mutation_counts = mtMut_counts_gene_cell,
      mutation_per_MB = Mut_MB_cell, mutation_per_copy = Mut_copy_cell, 
      MSS = MSS_gene_cell, MSS_oe = MSS_oe, MSS_weighted = MSS_gene_cell_weighted
    )
  }
  
  df2 <- do.call(rbind, lapply(1:37, function(i) {
    calculate_gene_metrics(i, mtGene$Starting[i]:mtGene$Ending[i], mtGene$Symbol[i])
  }))
  
  df3 <- do.call(rbind, lapply(1:6, function(i) {
    calculate_gene_metrics(i, unlist(mtRegion[i]), names(mtRegion)[i])
  }))
  
  df <- rbind(df1, df2, df3)
  return(df)
}


# quantify_MMB_psedobulk

quantify_MMB_psedobulk <- function(mmat, cov_per_pos_cell, SampleID){
  var_meta <- left_join(as.data.frame(rowData(mmat)), MLC_score, by = c("position" = "Position"))
  
  # Calculate per position counts and MLC score
  mtCov_pos_psedobulk <- rowSums(cov_per_pos_cell)
  mtMut_counts_psedobulk <- rowSums(assays(mmat)$allele_frequency * assays(mmat)$coverage)
  mtCov_psedobulk <- sum(mtCov_pos_psedobulk)
  mtCopy_psedobulk <- mtCov_psedobulk / 16569
  mtMut_psedobulk <- sum(mtMut_counts_psedobulk)
  Mut_MB_psedobulk <- mtMut_psedobulk / mtCov_psedobulk * 1e6
  Mut_copy_psedobulk <- mtMut_psedobulk / mtCopy_psedobulk
  MSS <- sum(var_meta$MLC_pos_score)
  MSS_weighted <- sum(var_meta$MLC_pos_score * var_meta$mean)
  
  df1 <- data.frame(barcode = "psedobulk",
                    sample = SampleID, symbol = "Total", 
                    coverage = mtCov_psedobulk, copy_number_abs = mtCopy_psedobulk, mutation_counts = mtMut_psedobulk,
                    mutation_per_MB = Mut_MB_psedobulk, mutation_per_copy = Mut_copy_psedobulk,
                    MSS = MSS, MSS_oe = MSS / sum(MLC_score$MLC_pos_score), MSS_weighted = MSS_weighted
  )
  
  calculate_metrics <- function(region, symbol) {
    row_idx <- var_meta$position %in% region
    mtCov_gene_psedobulk <- sum(mtCov_pos_psedobulk[region])
    mtCopy_gene_psedobulk <- mtCov_gene_psedobulk / length(region)
    
    if (sum(row_idx) > 0) {
      mtMut_counts_gene_psedobulk <- sum(mtMut_counts_psedobulk[row_idx])
      MSS_gene_psedobulk <- sum(var_meta$MLC_pos_score[row_idx])
      MSS_gene_psedobulk_weighted <- sum((var_meta$MLC_pos_score * var_meta$mean)[row_idx])
      MSS_oe <- MSS_gene_psedobulk / sum(MLC_score$MLC_pos_score[MLC_score$Position %in% region])
    } else {
      mtMut_counts_gene_psedobulk <- MSS_gene_psedobulk <- MSS_gene_psedobulk_weighted <- MSS_oe <- 0
    }
    
    Mut_MB_psedobulk <- mtMut_counts_gene_psedobulk / mtCov_gene_psedobulk * 1e6
    Mut_copy_psedobulk <- mtMut_counts_gene_psedobulk / mtCopy_gene_psedobulk
    
    data.frame(barcode = "pseudobulk",
               sample = SampleID, symbol = symbol, 
               coverage = mtCov_gene_psedobulk, copy_number_abs = mtCopy_gene_psedobulk, 
               mutation_counts = mtMut_counts_gene_psedobulk,
               mutation_per_MB = Mut_MB_psedobulk, mutation_per_copy = Mut_copy_psedobulk, 
               MSS = MSS_gene_psedobulk, MSS_oe = MSS_oe, MSS_weighted = MSS_gene_psedobulk_weighted
    )
  }
  
  df2 <- do.call(rbind, lapply(1:37, function(i) {
    calculate_metrics(mtGene$Starting[i]:mtGene$Ending[i], mtGene$Symbol[i])
  }))
  
  df3 <- do.call(rbind, lapply(1:6, function(i) {
    calculate_metrics(unlist(mtRegion[i]), names(mtRegion)[i])
  }))
  
  df <- rbind(df1, df2, df3)
  return(df)
}
