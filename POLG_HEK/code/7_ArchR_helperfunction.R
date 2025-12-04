cutoff_FDR  <- 0.05
cutoff_LFC  <- 0.6

up_labels <- c(
  # ECM / adhesion / mechanosensing
  "CYR61", "TGFBI", "ADAMTS1", 
  # cell–cell adhesion / protocadherins
  "PCDHA4", "PCDHA5", "PCDHA8", "PCDHA9",
  # KRAB-ZNF / transcriptional rewiring (just a few reps)
  "ZNF462", "ZNF433", "ZNF544", "ZNF813"
)

down_labels <- c(
  # mito genome maintenance & OXPHOS
  "MGME1", "MTO1", "PET117", "TIMMDC1", "COX8A", "TYMP", "SLC25A31", "SLC25A53"
)

highlight_genes  <- c(up_labels, down_labels)



plotVP <- function(df){
  df_volcano <- df %>%
    mutate(
      status = case_when(
        FDR <= cutoff_FDR & Log2FC >=  cutoff_LFC  ~ "up",
        FDR <= cutoff_FDR & Log2FC <= -cutoff_LFC  ~ "down",
        TRUE                                       ~ "ns"
      ),
      negLog10FDR = -log10(FDR)
    )
  
  ggplot(df_volcano, aes(x = Log2FC, y = negLog10FDR)) +
    geom_point(aes(color = status), alpha = 0.7, size = 1.8) +
    scale_color_manual(
      values = c("up" = "red", "down" = "blue", "ns" = "grey80"),
      name   = "Regulation",
      labels = c("down" = "Downregulated", "ns" = "Not significant", "up" = "Upregulated")
    ) +
    geom_vline(xintercept = c(-cutoff_LFC, cutoff_LFC), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(cutoff_FDR), linetype = "dashed", color = "grey60") +
    geom_text_repel(
      data = df_volcano %>% filter(status != "ns") %>% 
        filter(name %in% highlight_genes),              # <- only highlight genes
      # if you want only significant ones, use:
      # filter(name %in% highlight_genes, status != "ns"),
      aes(label = name, fontface = "bold.italic"),
      nudge_x = .15,
      box.padding = 0.5,
      nudge_y = 1,
      segment.curvature = -0.1,
      segment.ncp = 3,
      segment.angle = 20, size = 2.5, max.overlaps = 20
    ) +
    labs(
      x = "log2 fold change",
      y = "-log10(FDR)",
    ) +
    theme_classic() +
    theme(legend.position = "none")
}

