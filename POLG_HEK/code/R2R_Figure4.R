setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/deep_seq/")
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(SummarizedExperiment)
library(data.table)
'%ni%' <- Negate('%in%')
color_vec <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")

cov.list <- list(readRDS("output/2_mtDNA_pos_coverage.CTRL.rds"),
                 readRDS("output/2_mtDNA_pos_coverage.KI36.rds"),
                 readRDS("output/2_mtDNA_pos_coverage.KIA2.rds"))

cov.df <- as.data.frame(do.call(cbind, lapply(cov.list, rowMeans))) 
colnames(cov.df) <- c("CTRL", "KI36", "KIA2")
MLC <- read.table("../../database/MLC_score.tsv", header = TRUE)
# MITOMAP disease variants 
mitomap.disease <- read.table(
  "../../database/MITOMAP_disease_2022-05-25.txt",
  header = TRUE, sep = "\t", na.strings = c("", "NA"),
  fill = TRUE, comment.char = "", quote = ""
) %>%
  mutate(variant = paste0(pos, ref, ">", alt)) %>%
  filter(status %in% c("Cfrm", "Cfrm [P*]", " Cfrm [LP*]")) %>%  
  pull(pos) %>% unique()

cov.df <- cov.df %>% 
  mutate(Position = 1:n(), 
         pathogenic = ifelse(Position %in% mitomap.disease, T, F)) %>%
  left_join(MLC, by = "Position") %>%
  pivot_longer(cols = CTRL:KIA2, names_to = "condt", values_to = "cov")

cov.df %>% 
  ggplot(aes(x = condt, y = cov, fill = pathogenic)) +
  geom_violin(scale = "width") +
  theme_classic() +
  ylab("Mean Coverage") +
  xlab("Condition") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey"))

cov.df %>%
  mutate(
    MLC_bin = cut(
      MLC_pos_score,
      breaks = c(0, 0.33, 0.66, 1),
      labels = c("Low (0–0.33)", "Mid (0.33–0.66)", "High (0.66–1)"),
      include.lowest = TRUE
    )
  ) %>%
  ggplot(aes(x = condt, y = cov, fill = MLC_bin)) +
  geom_violin(scale = "width", trim = FALSE) +
  theme_classic() +
  ylab("Mean Coverage") +
  xlab("Condition") +
  scale_fill_manual(
    values = c(
      "Low (0–0.33)"  = "pink",
      "Mid (0.33–0.66)" = "red3",
      "High (0.66–1)" = "darkred"
    )
  ) +
  labs(fill = "MLC score interval")


