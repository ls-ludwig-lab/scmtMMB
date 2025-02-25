setwd("~/Ludwig_lab/scmtMMB/PBMC/code/")
library(ggplot2)
library(scales)
library(dplyr)
'%ni%' <- Negate('%in%')


scmtMMB <- read.csv("../PBMC_large_data_files/output/2_scmtMMB.csv", row.names = "X")
scmtMMB$X3243 <- "+"
scmtMMB_ex <- read.csv("../PBMC_large_data_files/output/2_scmtMMB_exclude_3243.csv", row.names = "X") 
scmtMMB_ex$X3243 <- "-"

color_vec <- c("M80" = "darkorchid4","M60" = "firebrick4", "M35" = "violetred3", 
               "M29" = "firebrick1", "H47" = "blue3", "H05" = "dodgerblue3")



# Define a function to handle the quantile filtering
apply_quantile_filter <- function(data, value_col, lower_quantile = 0.05, upper_quantile = 0.95) {
  data %>%
    group_by(sample, symbol) %>%
    mutate(
      q5 = quantile({{ value_col }}, lower_quantile),
      q95 = quantile({{ value_col }}, upper_quantile)
    ) %>%
    filter({{ value_col }} >= q5, {{ value_col }} <= q95) %>%
    ungroup()
}

# Define a function to create violin plots
create_violin_plot <- function(data, value_col, y_label, file_path, facet_ncol = 5, width = 10, height = 6, is_log = 2) {
  p <- data %>%
    ggplot(aes(x = sample, y = {{ value_col }})) +
    geom_violin(aes(fill = sample), scale = "width", alpha = 0.8) +
    scale_fill_manual(values = color_vec) +
    geom_point(data = data %>% filter(barcode == "pseudobulk"),
               aes(group = sample), shape = 23, size = 1, fill = "white") +
    facet_wrap(~symbol, ncol = facet_ncol) +
    theme_bw() +
    xlab(NULL) +
    theme_1_nogrid +
    labs(y = y_label) +
    Seurat::NoLegend()
  
  if (is_log == 2) {
    p <- p +
      scale_y_continuous(trans = log2_trans(),
                         breaks = trans_breaks("log2", function(x) 2^x),
                         labels = trans_format("log2", math_format(2^.x)))
  }else{
    p <- p +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                    labels = trans_format("log10", math_format(10^.x)))
  }
  
  ggsave(plot = p, file_path, width = width, height = height)
}

# Shared theme
theme_1_nogrid <- theme(
  aspect.ratio = 1 / 1,
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

# Weighted MSS Visualization
p112_data <- scmtMMB %>%
  filter(symbol %in% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA")) %>%
  apply_quantile_filter(MSS_weighted)
create_violin_plot(
  data = p112_data %>% mutate(MSS_weighted = MSS_weighted + 0.01),
  value_col = MSS_weighted,
  y_label = bquote(paste(log[2], "(scwMSS + 0.01)")),
  file_path = "../plot/24_scwMSS_complex_MELAS_q5_q95.pdf",
  facet_ncol = 6,
  width = 12,
  height = 2
)

# Gene-specific MSS Visualization
p13_data <- scmtMMB_ex %>%
  filter(symbol %ni% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA", "Total") & !grepl("MT_T", symbol)) %>%
  apply_quantile_filter(MSS_weighted)
create_violin_plot(
  data = p13_data %>% mutate(MSS_weighted = MSS_weighted + 0.01),
  value_col = MSS_weighted,
  y_label = bquote(paste(log[2], "(scwMSS + 0.01)")),
  file_path = "../plot/24_scwMSS_gene_MELAS_ex3243_q5_q95.pdf",
  facet_ncol = 5
)

# ND3 and ND4L MSS Visualization
p14_data <- scmtMMB_ex %>%
  filter(symbol %in% c("MT_ND3", "MT_ND4L")) %>%
  apply_quantile_filter(MSS_weighted)
create_violin_plot(
  data = p14_data %>% mutate(MSS_weighted = MSS_weighted + 0.01),
  value_col = MSS_weighted,
  y_label = bquote(paste(log[2], "(scwMSS + 0.01)")),
  file_path = "../plot/24_scwMSS_ND3_4L_MELAS_ex3243_q5_q95.pdf",
  facet_ncol = 2,
  width = 4,
  height = 2
)

# Mutation per MB Visualization
p212_data <- scmtMMB %>%
  filter(symbol %in% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA")) %>%
  apply_quantile_filter(mutation_per_MB)
create_violin_plot(
  data = p212_data %>% mutate(mutation_per_MB = mutation_per_MB + 1),
  value_col = mutation_per_MB,
  y_label = bquote(paste(log[10], "(mtscMPM + 1)")),
  file_path = "../plot/24_scMMB_complex_MELAS_q5_95.pdf",
  facet_ncol = 6,
  width = 12,
  height = 2,
  is_log = 10
)

p23_data <- scmtMMB_ex %>%
  filter(symbol %ni% c("tRNA", "CIII", "CV", "CIV", "CI", "rRNA", "Total") & !grepl("MT_T", symbol)) %>%
  apply_quantile_filter(mutation_per_MB)
create_violin_plot(
  data = p23_data %>% mutate(mutation_per_MB = mutation_per_MB + 1),
  value_col = mutation_per_MB,
  y_label = bquote(paste(log[10], "(mtscMPM + 1)")),
  file_path = "../plot/24_scMMB_gene_MELAS_ex3243_q5_q95.pdf",
  facet_ncol = 5,
  is_log = 10
)

p24_data <- scmtMMB_ex %>%
  filter(symbol %in% c("MT_ND4L", "MT_ND3")) %>%
  apply_quantile_filter(mutation_per_MB)
create_violin_plot(
  data = p24_data %>% mutate(mutation_per_MB = mutation_per_MB + 1),
  value_col = mutation_per_MB,
  y_label = bquote(paste(log[10], "(mtscMPM + 1)")),
  file_path = "../plot/24_scMMB_ND3_4L_MELAS_ex3243_q5_q95.pdf",
  facet_ncol = 2,
  width = 4,
  height = 2,
  is_log = 10
)
