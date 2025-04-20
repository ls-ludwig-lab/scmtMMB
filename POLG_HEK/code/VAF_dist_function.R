# ---------- 3. Functions - entropy ----------

prepare_mtVAR_meta <- function(mmat, mito.annotation) {
  as.data.frame(rowData(mmat)) %>%
    filter(vmr >= 0.01) %>%
    left_join(mito.annotation, by = "variant") %>%
    arrange(desc(mean)) %>%
    mutate(
      rank = 1:n(),
      label = ifelse(mean > 0.01 & effect %in% c("truncating", "missense"), variant, NA)
    )
}

extract_variant_matrices <- function(meta_df, assay_data) {
  get_matrix <- function(variant_list) {
    as.matrix(assay_data[rownames(assay_data) %in% variant_list, ])
  }
  
  list(
    synonymous  = get_matrix(meta_df %>% filter(effect == "synonymous") %>% pull(variant)),
    missense    = get_matrix(meta_df %>% filter(effect == "missense") %>% pull(variant)),
    truncating  = get_matrix(meta_df %>% filter(effect == "truncating") %>% pull(variant)),
    disease     = get_matrix(disease.var)
  )
}

plot_top_vaf_distribution <- function(vaf_matrix, title_label) {
  variant_var <- apply(vaf_matrix, 1, var)
  top_variants <- names(sort(variant_var, decreasing = TRUE))[1:100]
  
  vaf_long <- melt(vaf_matrix)
  colnames(vaf_long) <- c("variant", "Cell", "VAF")
  vaf_long %>% mutate(Position = as.integer(sub("([0-9]+).*", "\\1", variant))
  ) %>%
    left_join(MLC, by = "Position") 
  vaf_top <- subset(vaf_long, variant %in% top_variants)
  vaf_top$variant <- factor(vaf_top$variant, levels = top_variants)
  
  p_hist <- ggplot(vaf_top, aes(x = VAF, fill = variant)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, position = "identity", alpha = 0.3) +
    labs(title = paste("VAF Histogram (", title_label, ")"), x = "VAF", y = "Density") +
    theme_bw() + ylim(0, 20) + theme(legend.position = "none")
  
  p_ecdf <- ggplot(vaf_top, aes(x = VAF, group = variant)) +
    stat_ecdf(geom = "step") +
    labs(title = paste("VAF ECDF (", title_label, ")"), x = "VAF", y = "Cumulative Probability") +
    theme_bw() + theme(legend.position = "none")
  
  return(list(hist = p_hist, ecdf = p_ecdf))
}

calculate_entropy <- function(matrix, bins = 20) {
  apply(matrix, 1, function(vafs) {
    h <- hist(vafs, breaks = bins, plot = FALSE)
    entropy(h$counts, unit = "log2")
  })
}

plot_entropy_comparison <- function(group_matrices) {
  entropy_list <- lapply(group_matrices, calculate_entropy)
  boxplot(entropy_list,
          names = c("Synonymous", "Missense", "Truncating", "Disease"),
          main = "Shannon Entropy per Variant", ylab = "Entropy (bits)")
  print(kruskal.test(entropy_list))
}

plot_rank_ordered <- function(meta_df, title = "") {
  ggplot(meta_df, aes(x = rank, y = mean * 100, label = label, color = effect)) +
    geom_point() +
    ggrepel::geom_label_repel() +
    theme_classic() +
    scale_color_manual(values = c(
      "synonymous" = "grey50",
      "non-coding" = "grey50",
      "truncating" = "red",
      "missense" = "blue"
    )) +
    labs(x = "Rank ordered mutations", y = "Pseudobulk heteroplasmy (%)", title = title)
}

if (!dir.exists("../plot/VAF_dist_plots")) dir.create("../plot/VAF_dist_plots")

save_plot_pair <- function(plot_list, base_name) {
  ggsave(filename = paste0("../plot/VAF_dist_plots/", base_name, "_histogram.pdf"), plot = plot_list$hist, width = 2, height = 2)
  ggsave(filename = paste0("../plot/VAF_dist_plots/", base_name, "_ecdf.pdf"), plot = plot_list$ecdf, width = 2, height = 2)
}

save_rank_plot <- function(meta_df, sample_id) {
  p <- plot_rank_ordered(meta_df, paste0(sample_id, " - Rank ordered mutations"))
  ggsave(filename = paste0("../plot/VAF_dist_plots/", sample_id, "_rank_plot.pdf"), plot = p, width = 7, height = 5)
}

save_entropy_plot <- function(group_matrices, sample_id) {
  pdf(paste0("../plot/VAF_dist_plots/", sample_id, "_entropy_boxplot.pdf"), width = 5.5, height = 5.5)
  plot_entropy_comparison(group_matrices)
  dev.off()
}


# ---- Helper Functions ----
# Function to plot UMAP colored by a given variant level and save the plot as PDF
plot_variant <- function(df, variant, title, sampleName, umap_w = 4) {
  
  # Create UMAP plot
  p1 <- df %>%
    arrange(!!sym(variant)) %>%
    ggplot(aes(x = atacUMAP_1, y = atacUMAP_2, color = .data[[variant]])) +
    geom_point(size = 0.5) +
    scale_color_viridis_c(limits = quantile(df[[variant]], prob = c(0.05, 0.95))) +
    theme_void() +
    ggtitle(title)
  
  # Save UMAP plot as PDF
  pdf_file_umap <- paste0("../plot/VAF_dist_plots/", sampleName, "_", variant, "_UMAP.pdf")
  ggsave(filename = pdf_file_umap, plot = p1, device = "pdf", width = umap_w, height = 2)
  
  # Create Histogram plot
  p2 <- df %>%
    arrange(!!sym(variant)) %>%
    ggplot(aes(x = .data[[variant]] * 100)) +
    geom_histogram(binwidth = 10, fill = "white", color = "black") +
    theme_classic() +
    scale_y_continuous(limits = c(0, NA), expand = c(0,0)) +
    xlab(NULL) +
    ggtitle(title)
  
  # Save histogram plot as PDF
  pdf_file_hist <- paste0("../plot/VAF_dist_plots/", sampleName, "_", variant, "_hist.pdf")
  ggsave(filename = pdf_file_hist, plot = p2, device = "pdf", width = 4, height = 2)
  
  # Optionally print the plots to the active device
  print(p1)
  print(p2)
}



# Helper function to compute observed/expected co-occurrence with permutation p-value
test_independence_pair <- function(df_bin, var1, var2, n_perm = 1000) {
  N <- nrow(df_bin)
  A <- df_bin[[var1]]
  B <- df_bin[[var2]]
  
  obs_A <- sum(A)
  obs_B <- sum(B)
  obs_AB <- sum(A & B)
  
  pA <- obs_A / N
  pB <- obs_B / N
  expected_AB <- pA * pB * N
  log10_OE <- log10((obs_AB + 1e-10) / (expected_AB + 1e-10))  # Add small constant to avoid log2(0)
  
  # Permutation test (two-sided)
  perm_vals <- replicate(n_perm, {
    B_perm <- sample(B)
    sum(A & B_perm)
  })
  
  p_perm <- mean(abs(perm_vals - expected_AB) >= abs(obs_AB - expected_AB))
  
  return(data.frame(
    Variant1 = var1,
    Variant2 = var2,
    Observed = obs_AB,
    Expected = expected_AB,
    log10_OE = log10_OE,
    P_perm = p_perm,
    Significance = ifelse(p_perm < 0.05, "Significant", "Not Significant")
  ))
}


# Main function to run pairwise independence tests by variant class
run_independence_analysis_by_class <- function(se_object, meta_df, top_n = 100, n_perm = 1000, save_path = "./indep_results", dataset_label = "Dataset") {
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  
  variant_classes <- c("synonymous", "missense", "pathogenic", "truncating")
  results_all <- list()
  
  for (effect in variant_classes) {
    message("Processing: ", effect)
    variant_ids <- meta_df %>%
      filter(effect == !!effect) %>%
      arrange(desc(mean)) %>%
      slice_head(n = top_n) %>%
      pull(variant)
    
    if (length(variant_ids) < 2) {
      message("Not enough variants for ", effect)
      next
    }
    
    df_bin <- as.data.frame(t(assay(se_object)[variant_ids, ]) > 0)
    pairs <- combn(variant_ids, 2, simplify = FALSE)
    
    res <- purrr::map_dfr(pairs, ~test_independence_pair(df_bin, .x[1], .x[2], n_perm = n_perm)) %>%
      mutate(Variant_Type = effect)
    
    results_all[[effect]] <- res
    
    # Save CSV and heatmap for each type
    write.csv(res, file = file.path(save_path, paste0(dataset_label, "_", effect, "_indep_results.csv")), row.names = FALSE)
    
    heatmap_plot <- plot_independence_heatmap(res, paste("log10(O/E) Heatmap:", effect))
    ggsave(filename = file.path(save_path, paste0(dataset_label, "_", effect, "_heatmap.pdf")), plot = heatmap_plot, width = 6, height = 5)
  }
  
  # Save boxplot
  boxplot_all <- plot_independence_boxplot(results_all, title = paste("log10(O/E) by Variant Type -", dataset_label))
  ggsave(filename = file.path(save_path, paste0(dataset_label, "_log10OE_boxplot.pdf")), plot = boxplot_all, width = 5, height = 3)
  
  # Save bar plot of relationship types
  all_combined <- bind_rows(results_all, .id = "Variant_Type") %>%
    mutate(Relationship = case_when(
      log10_OE > 0.3 ~ "Co-occurring",
      log10_OE < -0.3 ~ "Exclusive",
      TRUE ~ "Neutral"
    ))
  
  bar_plot <- ggplot(all_combined, aes(x = Variant_Type, fill = Relationship)) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() +
    labs(title = paste("Proportion of Relationships by Variant Type -", dataset_label),
         y = "Proportion", x = NULL)
  
  ggsave(filename = file.path(save_path, paste0(dataset_label, "_relationship_barplot.pdf")), plot = bar_plot, width = 6, height = 5)
  
  # Chi-squared test
  chi_table <- table(all_combined$Variant_Type, all_combined$Relationship)
  chisq_result <- chisq.test(chi_table)
  saveRDS(chisq_result, file = file.path(save_path, paste0(dataset_label, "_chi_squared_result.rds")))
  
  return(results_all)
}

# Visualization: heatmap of log2(O/E)
plot_independence_heatmap <- function(results_df, title, log10_or_limits = c(-1.5, 1.5), p_cutoff = 0.05) {
  library(reshape2)
  # results_df <- results_df %>%
  #   mutate(log10_OE_masked = ifelse(P_perm < p_cutoff, log10_OE, NA))
  
  mat <- reshape2::dcast(results_df, Variant1 ~ Variant2, value.var = "log10_OE")
  rownames(mat) <- mat$Variant1
  mat <- mat[, -1, drop = FALSE]
  full_mat <- as.matrix(mat)
  full_mat[lower.tri(full_mat)] <- t(full_mat)[lower.tri(full_mat)]
  diag(full_mat) <- NA
  
  melted <- reshape2::melt(full_mat, na.rm = TRUE, value.name = "log10_OE")
  
  ggplot(melted, aes(Var1, Var2, fill = log10_OE)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", midpoint = 0,
      limits = log10_or_limits, oob = scales::squish,
      name = "log10(O/E)"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = NULL, y = NULL)
}

# Visualization: boxplot of log10(O/E) by variant type
plot_independence_boxplot <- function(results_all, title = "") {
  all_combined <- bind_rows(results_all, .id = "Variant_Type") %>% 
    mutate(Variant_Type = factor(Variant_Type, levels = c("synonymous", "missense", "truncating", "pathogenic"))) 
  
  ggplot(all_combined, aes(x = Variant_Type, y = log10_OE)) +
    geom_boxplot(outlier.alpha = 1, outlier.size = 0.5) +
    coord_cartesian(ylim = c(-12, 3)) +
    theme_bw() +
    labs(title = title, y = "log10(O/E)", x = NULL)
}


plot_variant_pair_scatter <- function(Name, se_object, var1, var2, title = "Variant Co-occurrence Scatterplot") {
  # Extract VAFs for the two variants
  df <- as.data.frame(t(assay(se_object)[c(var1, var2), ]))
  
  # Count total cells and non-zero counts for each variant
  n_total <- nrow(df)
  
  n_x <- sum(df[[var1]] > 0)
  pct_x <- round(100 * n_x / n_total, 1)
  
  n_y <- sum(df[[var2]] > 0)
  pct_y <- round(100 * n_y / n_total, 1)
  
  # Format axis labels
  x_label <- paste0("m.", var1, " (n > 0: ", n_x, ", ", pct_x, "%)")
  y_label <- paste0("m.", var2, " (n > 0: ", n_y, ", ", pct_y, "%)")
  
  # Plot
  ggplot(df, aes(x = .data[[var1]], y = .data[[var2]])) +
    geom_point() +
    labs(
      x = x_label,
      y = y_label,
      title = paste0(
        Name, " ", title, " –\n", 
        var1, " vs ", var2
      )
    ) +
    theme_bw()
}

# indep_results2$truncating %>% arrange(desc(abs(log10_OE)))
# plot_variant_pair_scatter("KI36", mmat_list[[2]], "6381G>A", "12767G>A") 
# #6381G>A 12767G>A        0 10.1607916 -11.006928  0.000     Significant   truncating



plot_variant_pair_grid <- function(Name, se_object, variant_vector, output_dir = "../plot") {
  df <- as.data.frame(t(assay(se_object)[variant_vector, ]))
  n <- length(variant_vector)
  colnames(df) <- variant_vector
  
  plot_list <- list()
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      var1 <- variant_vector[j]
      var2 <- variant_vector[i]
      
      if (i == j) {
        p <- ggplot() + theme_void()
      } else if (i > j) {
        # Lower triangle: scatterplot
        data_pair <- df[, c(var1, var2)]
        n_total <- nrow(data_pair)
        
        n_x <- sum(data_pair[[var1]] > 0)
        pct_x <- round(100 * n_x / n_total, 1)
        
        n_y <- sum(data_pair[[var2]] > 0)
        pct_y <- round(100 * n_y / n_total, 1)
        
        pct_both <- round(100 * sum(data_pair[[var1]] > 0 & data_pair[[var2]] > 0) / n_total, 1)
        
        x_label <- paste0("m.", var1, "\n(n > 0: ", n_x, ", ", pct_x, "%)")
        y_label <- paste0("m.", var2, "\n(n > 0: ", n_y, ", ", pct_y, "%)")
        title_text <- paste0("Double+ cells: ", pct_both, "%")
        
        p <- ggplot(data_pair, aes(x = .data[[var1]], y = .data[[var2]])) +
          geom_point(alpha = 0.4, size = 0.8) +
          labs(
            x = if (i == n) x_label else NULL,
            y = if (j == 1) y_label else NULL,
            title = title_text
          ) +
          theme_bw(base_size = 8) +
          theme(
            axis.title.x = element_text(size = 7),
            axis.title.y = element_text(size = 7),
            axis.text.x = if (i == n) element_text(size = 6) else element_blank(),
            axis.text.y = if (j == 1) element_text(size = 6) else element_blank(),
            plot.title = element_text(size = 7)
          )
      } else {
        # Upper triangle: empty
        p <- ggplot() + theme_void()
      }
      
      plot_list[[length(plot_list) + 1]] <- p
    }
  }
  
  # Save to PDF
  pdf_file <- file.path(output_dir, paste0(gsub(" ", "_", Name), "_variant_pair_grid.pdf"))
  pdf(pdf_file, width = n * 1.8, height = n * 1.8)
  gridExtra::grid.arrange(grobs = plot_list, ncol = n,
                          top = paste(Name, "- Variant Pair Grid"))
  dev.off()
  
  message("Saved plot to: ", pdf_file)
}


plot_variant_pair_grid_clean <- function(Name, se_object, variant_vector, output_dir = "../plot") {
  df <- as.data.frame(t(assay(se_object)[variant_vector, ]))
  n <- length(variant_vector)
  colnames(df) <- variant_vector
  
  plot_list <- list()
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      var1 <- variant_vector[j]
      var2 <- variant_vector[i]
      
      if (i == j || i < j) {
        # Empty or upper triangle
        p <- ggplot() + theme_void()
      } else {
        data_pair <- df[, c(var1, var2)]
        n_total <- nrow(data_pair)
        
        pct_both <- round(100 * sum(data_pair[[var1]] > 0 & data_pair[[var2]] > 0) / n_total, 1)
        title_text <- paste0("Double+ cells: ", pct_both, "%")
        
        p <- ggplot(data_pair, aes(x = .data[[var1]], y = .data[[var2]])) +
          geom_point(alpha = 0.4, size = 0.8) +
          labs(title = title_text) +
          coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
          theme_bw(base_size = 8) +
          theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(size = 7, hjust = 0.5),
            panel.grid = element_blank()
          )
      }
      
      plot_list[[length(plot_list) + 1]] <- p
    }
  }
  
  pdf_file <- file.path(output_dir, paste0(gsub(" ", "_", Name), "_scatter_grid_clean.pdf"))
  pdf(pdf_file, width = n * 1.5, height = n * 1.5)
  gridExtra::grid.arrange(grobs = plot_list, ncol = n)
  dev.off()
  
  message("Saved scatter grid to: ", pdf_file)
}