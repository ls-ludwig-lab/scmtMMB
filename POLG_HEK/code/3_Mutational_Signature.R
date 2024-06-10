##################################
# code source: caleblareau/mtscATACpaper_reproducibility/figure_paired_cd34_pbmc/code/10_mutational_signature.R
##################################
setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
library(data.table)
library(dplyr)
library(Seurat)
library(Signac)
library(BuenColors)
library(ggh4x)
'%ni%' <- Negate('%in%')
source("../../global_func/variant_calling.R")

# Construct trinucleotide context
refallele <- read.csv("../../database/refallele.csv")

# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}

# Process 3 digit signature based on letters
colnames(refallele) <- c("pos", "ref")

refallele$ref <- toupper(refallele$ref)
l <- as.character(refallele$ref)

# Gs happen to be at the first and last position
refallele$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))

# Remove Ns
refallele <- refallele[!grepl("N", refallele$three),]

# Make every possible mutation
refallele_long <- rbind(refallele,refallele, refallele,refallele)
refallele_long$alt <- rep(c("A", "C", "G", "T"), each = dim(refallele)[1])
refallele_long <- refallele_long[refallele_long$ref != refallele_long$alt,]

# add some meta data
refallele_long$variant <- paste0(as.character(refallele_long$pos), refallele_long$ref, ">", refallele_long$alt)
refallele_long$change <- paste0(refallele_long$ref, refallele_long$alt)
refallele_long$change_rc <- reverse_complement(paste0(refallele_long$ref, refallele_long$alt))

# A/G rich strand is "heavy" -- https://en.wikipedia.org/wiki/Heavy_strand
table(refallele_long$ref) # so the reference strand is light (more C/T)
refallele_long$strand <- ifelse(refallele_long$ref %in% c("C","T"), "L", "H") # why?

# Change to C/T as ref allele
refallele_long$rc3 <- reverse_complement(refallele_long$three)
refallele_long$three_plot <- ifelse(refallele_long$strand == "L", refallele_long$three, refallele_long$rc3)
refallele_long$group_change <- ifelse(refallele_long$strand == "L", refallele_long$change, refallele_long$change_rc)

# load mtDNA variants
CTRL.mmat <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.CTRL.rds")
KI36.mmat <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KI36.rds")
KIA2.mmat <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KIA2.rds")
mtVar.bulk.df <- read.csv("../output/2_mtVar.bulk.meta.csv", row.names = "X")

CTRL_variants <- rownames(CTRL.mmat)
KI36_variants <- rownames(KI36.mmat)
KIA2_variants <- rownames(KIA2.mmat)
# 
# Parental <- mtVar.bulk.df %>% filter(occurrence == "Parental") %>% pull(variant) %>% unique()
# KI36_private <- mtVar.bulk.df %>% filter(occurrence == "Private" & condition == "KI36") %>% pull(variant)
# KIA2_private <- mtVar.bulk.df %>% filter(occurrence == "Private" & condition == "KIA2") %>% pull(variant)
# POLG_variants <- mtVar.bulk.df %>% filter(occurrence == "POLG-sahred" & condition == "KIA2") %>% pull(variant) %>% unique()

# Annotate with called variants
refallele_long$KIA2 <- refallele_long$variant %in% KIA2_variants
refallele_long$KI36 <- refallele_long$variant %in% KI36_variants
refallele_long$CTRL <- refallele_long$variant %in% CTRL_variants
# refallele_long$Parental <- refallele_long$variant %in% Parental
# refallele_long$KI36P <- refallele_long$variant %in% KI36_private
# refallele_long$KIA2P <- refallele_long$variant %in% KIA2_private
# refallele_long$POLG <- refallele_long$variant %in% POLG_variants


total <- dim(refallele_long)[1]

total_KIA2 <- sum(refallele_long$KIA2)
total_CTRL <- sum(refallele_long$CTRL)
total_KI36 <- sum(refallele_long$KI36)
# total_Parental <- sum(refallele_long$Parental)
# total_KI36P <- sum(refallele_long$KI36P)
# total_KIA2P <- sum(refallele_long$KIA2P)
# total_POLG <- sum(refallele_long$POLG)

prop_df <- refallele_long %>% group_by(three_plot, group_change, strand) %>%
  summarize(observed_prop_KIA2 = sum(KIA2)/total_KIA2,
            observed_prop_KI36 = sum(KI36)/total_KI36,
            observed_prop_CTRL = sum(CTRL)/total_CTRL,
            # observed_prop_Parental = sum(Parental)/total_Parental,
            # observed_prop_KI36P = sum(KI36P)/total_KI36P,
            # observed_prop_KIA2P = sum(KIA2P)/total_KIA2P,
            # observed_prop_POLG = sum(POLG)/total_POLG,
            expected_prop = n()/total, n = n()) %>%
  mutate(fc_KIA2 = observed_prop_KIA2/expected_prop,
         fc_KI36 = observed_prop_KI36/expected_prop,
         fc_CTRL = observed_prop_CTRL/expected_prop,
         # fc_Parental = observed_prop_Parental/expected_prop,
         # fc_KI36P = observed_prop_KI36P/expected_prop,
         # fc_KIA2P = observed_prop_KIA2P/expected_prop,
         # fc_POLG = observed_prop_POLG/expected_prop
         )

prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)

prop_df$group_change_plot <- gsub('^([A-Z]{1})([A-Z]+)$', '\\1>\\2', prop_df$group_change)

write.table(refallele_long, "../output/3_mutational_signature_refallele_long.csv")
write.table(prop_df, "../output/3_mutational_signature_prop_df.csv")

# Visualize

CTRL.plot <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_CTRL, ylim(0,8))) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot() + L_border() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0), limits = c(0,7)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate (Observed/Expected)") +
  facet_nested(~group_change_plot, scales="free_x", space="free") +
  ggtitle("Mutational Signature CTRL") 

ggsave(plot = CTRL.plot, width = 10, height = 5, filename = "../plot/31_Mutational_Signature_CTRL.pdf")


KI36.plot <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_KI36, ylim(0,8))) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot() + L_border() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0), limits = c(0,7)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate (Observed/Expected)") +
  facet_nested(~group_change_plot, scales="free_x", space="free") +
  ggtitle("Mutational Signature KI36", subtitle = "Cosine Similarity (CTRL) = 0.814")

ggsave(plot = KI36.plot, width = 10, height = 5, filename = "../plot/31_Mutational_Signature_KI36.pdf")

KIA2.plot <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_KIA2, ylim(0,8))) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot() + L_border() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0), limits = c(0,7)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate (Observed/Expected)") +
  facet_nested(~group_change_plot, scales="free_x", space="free") +
  ggtitle("Mutational Signature KIA2", subtitle = "Cosine Similarity (CTRL) = 0.943")

ggsave(plot = KIA2.plot, width = 10, height = 5, filename = "../plot/31_Mutational_Signature_KIA2.pdf")

# Cosine similarity
lsa::cosine(prop_df$observed_prop_CTRL, prop_df$observed_prop_KI36)
lsa::cosine(prop_df$observed_prop_CTRL, prop_df$observed_prop_KIA2)
lsa::cosine(prop_df$observed_prop_KI36, prop_df$observed_prop_KIA2)





# lsa::cosine(prop_df$observed_prop_KI36, prop_df$observed_prop_KI36P)
# lsa::cosine(prop_df$observed_prop_KIA2, prop_df$observed_prop_KIA2P)
# ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_KI36P, ylim(0,8))) +
#   geom_bar(stat = "identity", position = "dodge") + pretty_plot() + L_border() + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())+
#   scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
#   theme(legend.position = "bottom") +
#   scale_y_continuous(expand = c(0,0), limits = c(0,7)) +
#   geom_hline(yintercept = 1, linetype =2, color = "black") +
#   labs(x = "Change in nucleotide", y = "Substitution Rate (Observed/Expected)") +
#   facet_nested(~group_change_plot, scales="free_x", space="free") +
#   ggtitle("Mutational Signature KI36P") 
# 
# ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_KIA2P, ylim(0,8))) +
#   geom_bar(stat = "identity", position = "dodge") + pretty_plot() + L_border() + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())+
#   scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
#   theme(legend.position = "bottom") +
#   scale_y_continuous(expand = c(0,0), limits = c(0,7)) +
#   geom_hline(yintercept = 1, linetype =2, color = "black") +
#   labs(x = "Change in nucleotide", y = "Substitution Rate (Observed/Expected)") +
#   facet_nested(~group_change_plot, scales="free_x", space="free") +
#   ggtitle("Mutational Signature KIA2P") 
# 
# ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_POLG, ylim(0,8))) +
#   geom_bar(stat = "identity", position = "dodge") + pretty_plot() + L_border() + 
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())+
#   scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
#   theme(legend.position = "bottom") +
#   scale_y_continuous(expand = c(0,0), limits = c(0,7)) +
#   geom_hline(yintercept = 1, linetype =2, color = "black") +
#   labs(x = "Change in nucleotide", y = "Substitution Rate (Observed/Expected)") +
#   facet_nested(~group_change_plot, scales="free_x", space="free") +
#   ggtitle("Mutational Signature POLG") 

