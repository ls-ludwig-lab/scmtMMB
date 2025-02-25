setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code")
library(data.table)
library(dplyr)
library(ggplot2)

mtGene <- fread("../../database/Coding_region.tsv") %>% mutate(Length = Ending - Starting + 1)
mtGene$Positions <- sapply(1:37, function(i) mtGene$Starting[i]:mtGene$Ending[i])
mtRegion <- sapply(unique(mtGene$Region), function(x) unique(unlist(mtGene$Positions[mtGene$Region == x])))
GenePos <- sapply(unique(mtGene$Symbol), function(x) unique(unlist(mtGene$Positions[mtGene$Symbol == x])))
MLC_score <- fread("../../database/MLC_score.tsv") # Mitochondrial local constraint (MLC) 

theme_nogrid <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

colours_com <- c("tRNA" = "magenta4", "rRNA" = "sienna4", 
                 "CIV" = "palegreen4", "CI" = "royalblue2", "CV" = "red3","CIII" = "tan",
                 "non-coding" = "grey50")

P1 <- mtGene %>% dplyr::select(Region, Length) %>% 
  group_by(Region) %>% 
  summarise(length = sum(Length)) %>%
  add_row(Region = "non-coding", length = 16569 - sum(.$length)) %>% 
  mutate(pct = length/16569 *100) %>% 
  ggplot(aes(x="", y=length, fill=Region)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values = colours_com) + 
  theme_void() +
  theme(
    aspect.ratio = 1 / 1,
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

P1.nl <- P1 + Seurat::NoLegend()
P1.leg <- ggpubr::get_legend(P1) %>% ggpubr::as_ggplot()

ggsave(plot = P1.nl, "../plot/mtDNA_comp.pdf", width=2, height=2)
ggsave(plot = P1.leg, "../plot/mtDNA_comp_leg.pdf", width=2, height=2)


colours <- c("tRNA" = "magenta4", 
             "MT-RNR1" = "sienna4", "MT-RNR2" = "sienna4", 
             "MT-ND1" = "palegreen4", "MT-ND2" = "palegreen4", "MT-CO1" = "royalblue2", 
             "MT-CO2" = "royalblue2", "MT-ATP8" = "red3", "MT-ATP6" = "red3", 
             "MT-CO3" = "royalblue2", "MT-ND3" = "palegreen4", "MT-ND4L" = "palegreen4", 
             "MT-ND4" = "palegreen4", "MT-ND5" = "palegreen4", "MT-ND6" = "palegreen4", "MT-CYB" = "tan",
             "Noncoding" = "grey")

P3 <- mtGene %>%
  dplyr::select(Symbol, Region, Length) %>%
  arrange(desc(Length)) %>%
  mutate(Symbol = sub("_", "-", Symbol)) %>% 
  filter(grepl("C", Region)) %>%
  mutate(Symbol = factor(Symbol, levels = Symbol)) %>%  # Adjust the order of the levels
  ggplot(aes(x = Symbol, y = Length, fill = Symbol)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = colours) +
  theme_classic() + theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0,0))

P3.nl <- P3 + Seurat::NoLegend()
P3.leg <- ggpubr::get_legend(P3) %>% ggpubr::as_ggplot()
ggsave(plot = P3.nl, "../plot/mtGene_length.pdf", width=4, height=4)

df <- data.frame(
  gene = rep(names(mtRegion), sapply(mtRegion, length)),
  Position = unlist(mtRegion)
)

df <- left_join(df, MLC_score,"Position")
df %>%filter(complete.cases(.)) %>% 
  group_by(gene) %>% 
summarise(
  Q1 = quantile(MLC_pos_score, 0.25),
  Q3 = quantile(MLC_pos_score, 0.75),
  median_val = median(MLC_pos_score)
)


P2 <- ggplot(df, aes(x = 1, MLC_pos_score)) +
  geom_violin(aes(fill = gene), scale = "width") + geom_boxplot(width=0.1) + 
  facet_wrap(~gene, ncol =6) + theme_bw() + theme_nogrid +
  scale_fill_manual(values = colours_com)+ Seurat::NoLegend()

P2.nl <- P2 + Seurat::NoLegend()
P2.leg <- ggpubr::get_legend(P2) %>% ggpubr::as_ggplot()
ggsave(plot = P2.nl, "../plot/MLC_score_distributon.pdf", width=4, height=2)


df2 <- data.frame(
  gene = rep(names(GenePos), sapply(GenePos, length)),
  Position = unlist(GenePos)
)

df2 <- left_join(df2, MLC_score,"Position") %>% mutate(gene = sub("_", "-", gene))


colours <- c("tRNA" = "magenta4", 
             "MT-RNR1" = "sienna4", "MT-RNR2" = "sienna4", 
             "MT-ND1" = "palegreen4", "MT-ND2" = "palegreen4", "MT-CO1" = "royalblue2", 
             "MT-CO2" = "royalblue2", "MT-ATP8" = "red3", "MT-ATP6" = "red3", 
             "MT-CO3" = "royalblue2", "MT-ND3" = "palegreen4", "MT-ND4L" = "palegreen4", 
             "MT-ND4" = "palegreen4", "MT-ND5" = "palegreen4", "MT-ND6" = "palegreen4", "MT-CYB" = "tan",
             "Noncoding" = "grey")

ggplot(df2 %>% filter(!grepl("MT-T", gene)), aes(x =1, y = MLC_pos_score)) +
  geom_violin(aes(fill = gene)) + geom_boxplot(width=0.1) + 
  facet_wrap(~gene, ncol =5) + theme_bw() + theme_nogrid +
  scale_fill_manual(values = colours)+ Seurat::NoLegend()

ggplot(df2 %>% filter(grepl("MT-T", gene)), aes(x =1, y = MLC_pos_score)) +
  geom_violin(fill = "magenta4") +
  facet_wrap(~gene, ncol =5) + theme_bw() + theme_nogrid + Seurat::NoLegend()

