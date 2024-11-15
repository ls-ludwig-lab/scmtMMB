library(data.table)
library(dplyr)
library(ggplot2)

mtGene <- fread("~/Ludwig_lab/scmtMMB/database/Coding_region.tsv") %>% mutate(Length = Ending - Starting + 1)
mtGene$Positions <- sapply(1:37, function(i) mtGene$Starting[i]:mtGene$Ending[i])
mtRegion <- sapply(unique(mtGene$Region), function(x) unique(unlist(mtGene$Positions[mtGene$Region == x])))
GenePos <- sapply(unique(mtGene$Symbol), function(x) unique(unlist(mtGene$Positions[mtGene$Symbol == x])))
MLC_score <- fread("~/Ludwig_lab/scmtMMB//database/MLC_score.tsv") # Mitochondrial local constraint (MLC) 

theme_nogrid <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

df <- data.frame(
  gene = rep(names(mtRegion), sapply(mtRegion, length)),
  Position = unlist(mtRegion)
)

df <- left_join(df, MLC_score,"Position")

colours_com <- c("tRNA" = "magenta4", "rRNA" = "sienna4", 
             "CIV" = "palegreen4", "CI" = "royalblue2", 
             "MT-CO2" = "royalblue2", "CV" = "red3","CIII" = "tan")

ggplot(df, aes(x = 1, MLC_pos_score)) +
  geom_violin(aes(fill = gene), scale = "width") + geom_boxplot(width=0.1) + 
  facet_wrap(~gene, ncol =6) + theme_bw() + theme_nogrid +
  scale_fill_manual(values = colours_com)+ Seurat::NoLegend()



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

