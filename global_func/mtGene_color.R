mitochondrial_data <- read.delim("../../database/Coding_region.tsv") %>% 
  mutate(Symbol = gsub("MT_", "", Symbol)) %>% 
  mutate(color = ifelse(Region == "tRNA", "tRNA", Symbol)) %>% 
  rbind(., data.frame(Symbol = rep("Non",4),
                      Starting = c(1, 5730, 8270, 16023),
                      Ending = c(576, 5760, 8294, 16569),
                      Region = rep("Noncoding",4),
                      color = rep("Noncoding",4)
  ))

mtGene_colors <- c("tRNA" = "magenta4", 
             "RNR1" = "mediumaquamarine", "RNR2" = "sienna4", 
             "ND1" = "magenta", "ND2" = "mediumblue", "CO1" = "olivedrab", 
             "CO2" = "orange2", "ATP8" = "orchid4", "ATP6" = "red3", 
             "CO3" = "royalblue2", "ND3" = "palegreen4", "ND4L" = "grey0", 
             "ND4" = "pink4", "ND5" = "yellow4", "ND6" = "steelblue4", "CYB" = "tan",
             "Noncoding" = "grey")
