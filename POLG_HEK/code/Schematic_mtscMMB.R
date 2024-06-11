library(dplyr)
library(ggplot2)
setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
MLC_score <- data.table::fread("../../database/MLC_score.tsv") # Mitochondrial local constraint (MLC) 

mitochondrial_data <- read.delim("../../database/Coding_region.tsv") %>% 
  mutate(Symbol = gsub("MT_", "", Symbol)) %>% 
  mutate(color = ifelse(Region == "tRNA", "tRNA", Symbol)) %>% 
  rbind(., data.frame(Symbol = rep("Non",4),
                      Starting = c(1, 5730, 8270, 16023),
                      Ending = c(576, 5760, 8294, 16569),
                      Region = rep("Noncoding",4),
                      color = rep("Noncoding",4)
  ))

colours <- c("tRNA" = "magenta4", 
             "RNR1" = "mediumaquamarine", "RNR2" = "sienna4", 
             "ND1" = "magenta", "ND2" = "mediumblue", "CO1" = "olivedrab", 
             "CO2" = "orange2", "ATP8" = "orchid4", "ATP6" = "red3", 
             "CO3" = "royalblue2", "ND3" = "palegreen4", "ND4L" = "grey0", 
             "ND4" = "pink4", "ND5" = "yellow4", "ND6" = "steelblue4", "CYB" = "tan",
             "Noncoding" = "grey")

# visualize coverage
(p1 <- MLC_score %>% ggplot() +
  geom_line(aes(x = Position, y = 1, color = MLC_pos_score), linewidth = 15) +
  scale_color_viridis_c(option = "magma") +
  geom_rect(data = mitochondrial_data, 
            aes(xmin = Starting, xmax = Ending, ymin = 0.8, ymax =0.85 , fill = color), alpha = 1) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  scale_fill_manual(values = colours) +
  guides(fill=guide_legend(ncol=2)) +
  coord_polar() + theme_void())

p2 <- p1 + Seurat::NoLegend()
leg <- ggpubr::get_legend(p1)
p3 <- ggpubr::as_ggplot(leg)
ggsave("../plot/Scheme_MLC.pdf", p2, width = 5, height = 5)
ggsave("../plot/Scheme_MLC_leg.pdf", p3, width = 5, height = 5)

set.seed(1234)
library(VGAM)
# Set the parameters for left-skewed distribution
shape <- 2  # Shape parameter for left skewness
scale <- 0.01  # Scale parameter for mean

# Generate the left-skewed sequence
left_skewed_sequence <- abs(rlgamma(n = 16568, shape = shape, scale = scale))
left_skewed_sequence[sample(1:16568, 16000)] <- NA
# left_skewed_sequence[sample(1:16568, 10)] <- 1

p4 <- ggplot() +
  geom_point(aes(x = 1:16568, y = left_skewed_sequence, color = left_skewed_sequence)) +
  scale_color_viridis_c(begin = 0, na.value = "white") +
  geom_rect(data = mitochondrial_data, 
            aes(xmin = Starting, xmax = Ending, ymin = -0.005, ymax =-0 , fill = color), alpha = 1) +
  scale_y_log10() +
  scale_y_continuous(limits = c(-0.05, NA), expand = c(0,0)) +
  scale_fill_manual(values = colours) +
  guides(fill=guide_legend(ncol=2)) +
  coord_polar() + theme_minimal()
  
p5 <- p4 + Seurat::NoLegend()
leg <- ggpubr::get_legend(p4)
p6 <- ggpubr::as_ggplot(leg)
ggsave("../plot/Scheme_heteroplasmy.pdf", p5, width = 5, height = 5)
ggsave("../plot/Scheme_heteroplasmy_leg.pdf", p6, width = 5, height = 5)

