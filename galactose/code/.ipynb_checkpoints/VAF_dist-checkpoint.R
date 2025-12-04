# =======================================
# Mitochondrial Variant Analysis Pipeline
# =======================================

### Part 1: Vaf_distribution

# ---------- 1. Load Libraries ----------
setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(entropy)
library(GGally)
library(plotly)
library(gridExtra)
library(SummarizedExperiment)
source("VAF_dist_function.R")

# ---------- 2. Load Data ----------
# Define sample names
samples <- c("CTRL", "KI36", "KIA2")

# Read mmat and bulk data for each sample
mmat_list <- setNames(lapply(samples, function(s) {
  readRDS(sprintf("../POLG_HEK_large_data_files/output/2_mmtx.%s.rds", s))
}), samples)

# MITOMAP disease variants
mitomap.disease <- read.table(
  "../../database/MITOMAP_disease_2022-05-25.txt",
  header = TRUE, sep = "\t", na.strings = c("", "NA"),
  fill = TRUE, comment.char = "", quote = ""
)
# MLC score
MLC <- read.table("../../database/MLC_score.tsv", header = TRUE)

# mitoVAR annotation
mito.annotation <- read.table("../../database/mitovar_annotation.txt", header = TRUE) %>%
  mutate(
    variant = gsub("MT_(\\d+)_([A-Z])/([A-Z])", "\\1\\2>\\3", Uploaded_variation),
    Position = as.integer(sub("([0-9]+).*", "\\1", variant))
  ) %>%
  left_join(MLC, by = "Position") %>%
  select(variant, Consequence, SYMBOL, BIOTYPE, Protein_position, Amino_acids, MLC_pos_score, SIFT, PolyPhen) %>%
  mutate(
    effect = case_when(
      Consequence %in% c("stop_retained_variant", "synonymous_variant") ~ "synonymous",
      Consequence %in% c("missense_variant") ~ "missense",
      Consequence %in% c("intergenic_variant", "non_coding_transcript_exon_variant") ~ "non-coding",
      TRUE ~ "truncating"),
    SIFT_score = ifelse(
      SIFT == "-",
      NA,
      as.numeric(sub(".*\\(([^)]+)\\)", "\\1", SIFT))),
    PolyPhen_score = ifelse(
      PolyPhen == "-",
      NA,
      as.numeric(sub(".*\\(([^)]+)\\)", "\\1", PolyPhen))),
  )

# Define confirmed pathogenic variants
disease.var <- mitomap.disease %>%
  mutate(variant = paste0(pos, ref, ">", alt)) %>%
  filter(status %in% c("Cfrm", "Cfrm [P*]", " Cfrm [LP*]")) %>%
  pull(variant)

# Tag pathogenic
mito.annotation <- mito.annotation %>%
  mutate(effect = ifelse(variant %in% disease.var, "pathogenic", effect))


# ---------- 4. Apply to Both Datasets ----------

# For mmat_list[[2]]
meta2 <- prepare_mtVAR_meta(mmat_list[[2]], mito.annotation)
matrix2 <- extract_variant_matrices(meta2, assay(mmat_list[[2]]))

# For mmat_list[[3]]
meta3 <- prepare_mtVAR_meta(mmat_list[[3]], mito.annotation)
matrix3 <- extract_variant_matrices(meta3, assay(mmat_list[[3]]))

# ---------- 5. Plotting ----------
# --------- Sample: KI36 (mmat_list[[2]]) ---------
plots2_disease <- plot_top_vaf_distribution(matrix2$disease, "KI36 - Disease")
plots2_trunc   <- plot_top_vaf_distribution(matrix2$truncating, "KI36 - Truncating")
plots2_miss    <- plot_top_vaf_distribution(matrix2$missense, "KI36 - Missense")
plots2_syn     <- plot_top_vaf_distribution(matrix2$synonymous, "KI36 - Synonymous")

save_plot_pair(plots2_disease, "KI36_disease")
save_plot_pair(plots2_trunc,   "KI36_truncating")
save_plot_pair(plots2_miss,    "KI36_missense")
save_plot_pair(plots2_syn,     "KI36_synonymous")

save_rank_plot(meta2, "KI36")
save_entropy_plot(matrix2[c("synonymous", "missense", "truncating", "disease")], "KI36")

# --------- Sample: KIA2 (mmat_list[[3]]) ---------
plots3_disease <- plot_top_vaf_distribution(matrix3$disease, "KIA2 - Disease")
plots3_trunc   <- plot_top_vaf_distribution(matrix3$truncating, "KIA2 - Truncating")
plots3_miss    <- plot_top_vaf_distribution(matrix3$missense, "KIA2 - Missense")
plots3_syn     <- plot_top_vaf_distribution(matrix3$synonymous, "KIA2 - Synonymous")

save_plot_pair(plots3_disease, "KIA2_disease")
save_plot_pair(plots3_trunc,   "KIA2_truncating")
save_plot_pair(plots3_miss,    "KIA2_missense")
save_plot_pair(plots3_syn,     "KIA2_synonymous")

save_rank_plot(meta3, "KIA2")
save_entropy_plot(matrix3[c("synonymous", "missense", "truncating", "disease")], "KIA2")

#### plot some UMAP
# Load Seurat object and add UMAP metadata
SO <- readRDS("../POLG_HEK_large_data_files/output/1_POLG_HEK_preprocesssed_seurat.rds")
SO <- AddMetaData(SO, metadata = as.data.frame(SO@reductions$atac.umap@cell.embeddings))

# ---- Data Loading & Preparation ----
# Define sample names
samples <- c("CTRL", "KI36", "KIA2")

# Create single-cell data frames from the rowData and add a 'HetMut' column
sc_df <- lapply(mmat_list, function(mmat) {
  as.data.frame(rowData(mmat))
})

# ---- Variant UMAP Plots ----

# -- For KI36 --
# Subset KI36 cells (using colnames from the mmat object) and apply a UMAP filter
KI36 <- SO[, colnames(mmat_list[["KI36"]])]
KI36 <- subset(KI36, atacUMAP_1 < 4 & atacUMAP_2 > 1.9)

# Define KI36 variants and their display names
variants_KI36 <- c("X15710C_T", "X13414G_A", "X10197G_A", "X8313G_A")
variant_names_KI36 <- c("m.15710C>T", "m.13414G>A", "m.10197G>A", "m.8313G>A")

# For each variant, add the corresponding assay values to the Seurat object.
# (Here we assume that the row names in the assay are like "15710C>T" so we translate "X15710C_T" accordingly.)
for(v in variants_KI36) {
  row_name <- sub("_", ">", sub("^X", "", v))
  KI36[[v]] <- assay(mmat_list[["KI36"]])[row_name, colnames(KI36)]
}
# Update the metadata data frame after adding new columns
KI36.df <- as.data.frame(KI36@meta.data)
# Create UMAP plots for each variant and save as PDFs
for(i in seq_along(variants_KI36)) {
  plot_variant(KI36.df, variants_KI36[i], variant_names_KI36[i], sampleName = "KI36")
}

# -- For KIA2 --
KIA2 <- SO[, colnames(mmat_list[["KIA2"]])]
KIA2 <- subset(KIA2, atacUMAP_1 < 4 & atacUMAP_2 < -0.5)
variants_KIA2 <- c("X5348C_A", "X583G_A")
variant_names_KIA2 <- c("m.5348C>A", "m.583G>A")

for(v in variants_KIA2) {
  row_name <- sub("_", ">", sub("^X", "", v))
  KIA2[[v]] <- assay(mmat_list[["KIA2"]])[row_name, colnames(KIA2)]
}

# Plot using Seurat's FeaturePlot and then using ggplot for consistency
KIA2.df <- as.data.frame(KIA2@meta.data)
for(i in seq_along(variants_KIA2)) {
  plot_variant(KIA2.df, variants_KIA2[i], variant_names_KIA2[i], sampleName = "KIA2",  umap_w = 3)
}


## Part 2: Independence analysis

indep_results2 <- run_independence_analysis_by_class(
  se_object = mmat_list[[2]],
  meta_df = meta2,
  top_n = 100,
  n_perm = 1000,
  save_path = "../output/indep/KI36",
  dataset_label = "KI36"
)

indep_results3 <- run_independence_analysis_by_class(
  se_object = mmat_list[[3]],
  meta_df = meta3,
  top_n = 100,
  n_perm = 1000,
  save_path = "../output/indep/KIA2",
  dataset_label = "KIA2"
)

plot_independence_boxplot(indep_results2)
plot_independence_boxplot(indep_results3)


#### just plot some scatter plot 

KIA2.Top5.p <- meta3 %>% filter(effect == "pathogenic") %>% arrange(desc(mean)) %>% slice_head(n = 5) %>% pull(variant)
KIA2.Top5.t <- meta3 %>% filter(effect == "truncating") %>% arrange(desc(mean)) %>% slice_head(n = 5) %>% pull(variant)

KI36.Top5.p <- meta2 %>% filter(effect == "pathogenic") %>% arrange(desc(mean)) %>% slice_head(n = 5) %>% pull(variant)
KI36.Top5.t <- meta2 %>% filter(effect == "truncating") %>% arrange(desc(mean)) %>% slice_head(n = 5) %>% pull(variant)


plot_variant_pair_grid(
  Name = "KIA2_Top5_Pathogenic",
  se_object = mmat_list[[3]],
  variant_vector = KIA2.Top5.p
)

plot_variant_pair_grid(
  Name = "KIA2_Top5_Truncating",
  se_object = mmat_list[[3]],
  variant_vector = KIA2.Top5.t
)

plot_variant_pair_grid(
  Name = "KI36_Top5_Pathogenic",
  se_object = mmat_list[[2]],
  variant_vector = KI36.Top5.p
)

plot_variant_pair_grid(
  Name = "KI36_Top5_Truncating",
  se_object = mmat_list[[2]],
  variant_vector = KI36.Top5.t
)

plot_variant_pair_grid_clean(
  Name = "KIA2_Top5_Pathogenic",
  se_object = mmat_list[[3]],
  variant_vector = KIA2.Top5.p
)

plot_variant_pair_grid_clean(
  Name = "KIA2_Top5_Truncating",
  se_object = mmat_list[[3]],
  variant_vector = KIA2.Top5.t
)

plot_variant_pair_grid_clean(
  Name = "KI36_Top5_Pathogenic",
  se_object = mmat_list[[2]],
  variant_vector = KI36.Top5.p
)

plot_variant_pair_grid_clean(
  Name = "KI36_Top5_Truncating",
  se_object = mmat_list[[2]],
  variant_vector = KI36.Top5.t
)


lapply(matrix2, function(mtx){
  hist(colSums(mtx >0)/nrow(mtx))
})


plot_variant_presence_histograms <- function(matrix_list, sample_name = "Sample") {
  names(matrix_list) <- names(matrix_list)  # ensure names are preserved
  
  plots <- lapply(names(matrix_list), function(var_type) {
    mtx <- matrix_list[[var_type]]
    df <- data.frame(
      CellID = colnames(mtx),
      Prop_Cells = colSums(mtx > 0) / nrow(mtx)
    )
    
    ggplot(df, aes(x = Prop_Cells)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "black") +
      xlim(0,1)+
      theme_minimal() +
      labs(
        title = paste0(sample_name, " – ", var_type, " variants"),
        x = "Proportion of Detected Variants per cell",
        y = "Cell counts"
      )
  })
  
  names(plots) <- names(matrix_list)
  return(plots)
}

plots_KI36 <- plot_variant_presence_histograms(matrix2, sample_name = "KI36")

plots_KIA2 <- plot_variant_presence_histograms(matrix3, sample_name = "KIA2")



plot_variant_presence <- function(matrix_list, sample_name = "Sample") {
  names(matrix_list) <- names(matrix_list)  # ensure names are preserved
  
  # Define custom x-axis limits by index
  custom_xlim <- list(
    c(0, 400),
    c(0, 500),
    c(0, 70),
    c(0, 8)
  )
  
  plots <- lapply(seq_along(matrix_list), function(i) {
    var_type <- names(matrix_list)[i]
    mtx <- matrix_list[[i]]
    
    df <- data.frame(
      CellID = colnames(mtx),
      counts = colSums(mtx > 0)
    )
    
    ggplot(df, aes(x = counts)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "black") +
      theme_minimal() +
      labs(
        title = paste0(sample_name, " – ", var_type, " variants"),
        x = "Detected Variants per Cell",
        y = "Cell Count"
      ) +
      xlim(custom_xlim[[i]])
  })
  
  names(plots) <- names(matrix_list)
  return(plots)
}
plots_KI36_c <- plot_variant_presence(matrix2, sample_name = "KI36")
plots_KIA2_c <- plot_variant_presence(matrix3, sample_name = "KIA2")
