setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code")
if(!require(ineq)) install.packages("ineq", dependencies=TRUE)
library(ineq)
library(Matrix)
library(ggplot2)
library(Seurat)
source("../../global_func/mtGene_color.R")


summary_df <- lapply(c("CTRL","KI36", "KIA2"), function(condt){
  set.seed(1234)
  SO <- readRDS(paste0("../output/scVAF/VAF_SO_",condt,".rds"))
  meta <- read.csv(paste0("../output/scVAF/",condt,".csv"), row.names = "X")
  mtx <- SO@assays$VAF@counts
  
  # Define the target number of rows
  target_rows <- 16569
  
  if (nrow(mtx) < target_rows) {
    # Calculate how many rows of zeros need to be added
    rows_to_add <- target_rows - nrow(mtx)
    
    # Create a matrix of zeros with the same number of columns
    zero_matrix <- matrix(0, nrow = rows_to_add, ncol = ncol(mtx))
    
    # Bind the original matrix with the zero matrix
    mtx <- rbind(mtx, zero_matrix)
  }
  
  # Define a function to calculate the mean and Gini index for each column
  calculate_stats <- function(column) {
    mean_value <- mean(column)
    sum_value <- sum(column)
    gini_value <- Gini(column, )
    return(c(mean = mean_value, sum = sum_value, gini = gini_value))
  }
  
  summary_df <- as.data.frame(do.call(rbind, lapply(as.data.frame(mtx), calculate_stats)))
  summary_df$clone <- as.factor(meta[rownames(summary_df),"seurat_clusters"])
  summary_df$dataset = condt

  summary_df
})

summary_df <- do.call(rbind, summary_df)

p <- ggplot(summary_df, aes(x = sum, y = gini, col = dataset)) + geom_point() +
  # ggplot(summary_df, aes(x = mean, y = gini, col = dataset)) + geom_point() +
  #scale_x_continuous(limits = c(0, 0.005)) +
  scale_y_continuous(limits = c(0.95, 1)) +
  scale_color_manual(values = c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090"))+
  theme_classic() + NoLegend()

ggsave(plot = p, filename = paste0("../plot/scVAF/gini.pdf"), width = 5, height = 5)


summary_df2 <- lapply(c("CTRL","KI36", "KIA2"), function(condt){
  set.seed(1234)
  SO <- readRDS(paste0("../output/scVAF/VAF_SO_",condt,".rds"))
  meta <- read.csv(paste0("../output/scVAF/",condt,".csv"), row.names = "X")
  mtx <- SO@assays$VAF@counts
  
  calculate_stats <- function(column) {
    mean_value <- mean(column)
    sum_value <- sum(column)
    gini_value <- Gini(column)
    return(c(mean = mean_value, sum = sum_value, gini = gini_value))
  }

  summary_df <- as.data.frame(do.call(rbind, lapply(as.data.frame(mtx), calculate_stats)))
  summary_df$clone <- as.factor(meta[rownames(summary_df),"seurat_clusters"])
  summary_df$dataset = condt

  summary_df
})

summary_df2 <- do.call(rbind, summary_df2)

p2 <- ggplot(summary_df2, aes(x = sum, y = gini, col = dataset)) + geom_point() +
  # ggplot(summary_df2, aes(x = mean, y = gini, col = dataset)) + geom_point() +
  #scale_x_continuous(limits = c(0, 0.005)) +
  #scale_y_continuous(limits = c(0.95, 1)) +
  scale_color_manual(values = c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090"))+
  theme_classic() + NoLegend()

ggsave(plot = p2, filename = paste0("../plot/scVAF/gini_2.pdf"), width = 5, height = 5)
