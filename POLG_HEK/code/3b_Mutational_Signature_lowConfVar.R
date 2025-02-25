##################################
# code source: caleblareau/mtscATACpaper_reproducibility/figure_paired_cd34_pbmc/code/10_mutational_signature.R
##################################
setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
library(data.table)
library(dplyr)
library(BuenColors)
library(ggh4x)
'%ni%' <- Negate('%in%')
source("../../global_func/variant_calling.R")

mitoSE.CTRL <- readRDS("../POLG_HEK_large_data_files/output/2_mgatk.CTRL.rds")
mitoSE.KI36 <- readRDS("../POLG_HEK_large_data_files/output/2_mgatk.KI36.rds")
mitoSE.KIA2 <- readRDS("../POLG_HEK_large_data_files/output/2_mgatk.KIA2.rds")

mmat.CTRL <- call_mutations_mgatk(mitoSE.CTRL)
mmat.KI36 <- call_mutations_mgatk(mitoSE.KI36)
mmat.KIA2 <- call_mutations_mgatk(mitoSE.KIA2)

# filter low confident variant
mtMut <- function(mmat, str_low, str_high, n_cells = 1, vmr.thres = 0.01){
  var <- as.data.frame(rowData(mmat))
  var <- var %>% subset(n_cells_conf_detected >= n_cells &
                          strand_correlation > str_low &
                          strand_correlation <= str_high &
                          vmr >= vmr.thres &
                          variant %ni% c("301A>C", "302A>C", "310T>C", "316G>C")) %>% pull(variant)
  mmat <- mmat[var,]
  return(mmat)
}

# analyze three cell lines separately and visualize selected variants
mmat.KI36.high <- mtMut(mmat.KI36, 0.65, 1)
rowData(mmat.KI36.high)$condition <- "KI36"
mmat.KIA2.high <- mtMut(mmat.KIA2, 0.65, 1)
rowData(mmat.KIA2.high)$condition <- "KIA2"

mmat.KI36.med <- mtMut(mmat.KI36, 0.45, 0.65)
rowData(mmat.KI36.med)$condition <- "KI36"
mmat.KIA2.med <- mtMut(mmat.KIA2, 0.45, 0.65)
rowData(mmat.KI36.med)$condition <- "KIA2"

mmat.KI36.neg <- mtMut(mmat.KI36, -1.1, 0)
rowData(mmat.KI36.neg)$condition <- "KI36"
mmat.KIA2.neg <- mtMut(mmat.KIA2, -1.1, 0)
rowData(mmat.KIA2.neg)$condition <- "KIA2"


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
KI36_variants.h <- rownames(mmat.KI36.high)
KIA2_variants.h <- rownames(mmat.KIA2.high)
KI36_variants.m <- rownames(mmat.KI36.med)
KIA2_variants.m <- rownames(mmat.KIA2.med)
KI36_variants.n <- rownames(mmat.KI36.neg)
KIA2_variants.n <- rownames(mmat.KIA2.neg)

# Annotate with called variants
refallele_long$KIA2.h <- refallele_long$variant %in% KI36_variants.h
refallele_long$KI36.h <- refallele_long$variant %in% KIA2_variants.h
refallele_long$KIA2.m <- refallele_long$variant %in% KI36_variants.m
refallele_long$KI36.m <- refallele_long$variant %in% KIA2_variants.m
refallele_long$KIA2.n <- refallele_long$variant %in% KI36_variants.n
refallele_long$KI36.n <- refallele_long$variant %in% KIA2_variants.n

total <- dim(refallele_long)[1]
total_KIA2.h <- sum(refallele_long$KIA2.h)
total_KI36.h <- sum(refallele_long$KI36.h)
total_KIA2.m <- sum(refallele_long$KIA2.m)
total_KI36.m <- sum(refallele_long$KI36.m)
total_KIA2.n <- sum(refallele_long$KIA2.n)
total_KI36.n <- sum(refallele_long$KI36.n)

prop_df <- refallele_long %>% group_by(three_plot, group_change, strand) %>%
  summarize(observed_prop_KIA2.h = sum(KIA2.h)/total_KIA2.h,
            observed_prop_KI36.h = sum(KI36.h)/total_KI36.h,
            observed_prop_KIA2.m = sum(KIA2.m)/total_KIA2.m,
            observed_prop_KI36.m = sum(KI36.m)/total_KI36.m,
            observed_prop_KIA2.n = sum(KIA2.n)/total_KIA2.n,
            observed_prop_KI36.n = sum(KI36.n)/total_KI36.n,
            expected_prop = n()/total, n = n()) %>%
  mutate(fc_KIA2.h = observed_prop_KIA2.h/expected_prop,
         fc_KI36.h = observed_prop_KI36.h/expected_prop,
         fc_KIA2.m = observed_prop_KIA2.m/expected_prop,
         fc_KI36.m = observed_prop_KI36.m/expected_prop,
         fc_KIA2.n = observed_prop_KIA2.n/expected_prop,
         fc_KI36.n = observed_prop_KI36.n/expected_prop)
prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)

prop_df$group_change_plot <- gsub('^([A-Z]{1})([A-Z]+)$', '\\1>\\2', prop_df$group_change)

#write.table(refallele_long, "../output/3_mutational_signature_refallele_long.csv")
#write.table(prop_df, "../output/3_mutational_signature_prop_df.csv")

# Visualize

KI36.plot <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_KI36.m, ylim(0,8))) +
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
  ggtitle("Mutational Signature KI36 with str-cor 0.45-0.65")#, subtitle = "Cosine Similarity (KI36.h) = 0.829")

ggsave(plot = KI36.plot, width = 4, height = 2.5, filename = "../plot/31_Mutational_Signature_KI36_lcv_45_65.pdf")

KIA2.plot <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_KIA2.m, ylim(0,8))) +
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
  ggtitle("Mutational Signature KIA2 with str-cor 0.45-0.65")#, subtitle = "Cosine Similarity (KIA2.h) = 0.953")

ggsave(plot = KIA2.plot, width = 4, height = 2.5, filename = "../plot/31_Mutational_Signature_KIA2_lcv_45_65.pdf")


KI36.plot.n <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_KI36.n, ylim(0,8))) +
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
  ggtitle("Mutational Signature KI36 with str-cor <0")#, subtitle = "Cosine Similarity (KI36.h) = 0.829")

ggsave(plot = KI36.plot.n, width = 4, height = 2.5, filename = "../plot/31_Mutational_Signature_KI36_lcv_neg.pdf")

KIA2.plot.n <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_KIA2.n, ylim(0,8))) +
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
  ggtitle("Mutational Signature KIA2 with str-cor <0")#, subtitle = "Cosine Similarity (KIA2.h) = 0.953")

ggsave(plot = KIA2.plot.n, width = 4, height = 2.5, filename = "../plot/31_Mutational_Signature_KIA2_lcv_neg.pdf")


# Cosine similarity
lsa::cosine(prop_df$observed_prop_KI36.h, prop_df$observed_prop_KI36.m) # 0.9667426
lsa::cosine(prop_df$observed_prop_KI36.h, prop_df$observed_prop_KI36.n) # 0.1910642
lsa::cosine(prop_df$observed_prop_KIA2.h, prop_df$observed_prop_KIA2.m) # 0.9551243
lsa::cosine(prop_df$observed_prop_KIA2.h, prop_df$observed_prop_KIA2.n) # 0.5564896
