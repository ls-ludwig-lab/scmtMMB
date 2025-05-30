setwd("~/Ludwig_lab/scmtMMB/POLG_HEK/code/")
library(Matrix)
library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(viridis)
library(scales)
'%ni%' <- Negate('%in%')
color_vec <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")


# Load the seurat object
mtVar.bulk.df <- read.csv("../output/2_mtVar.bulk.meta.csv", row.names = "X")
CTRL_dnV <- mtVar.bulk.df %>% filter(occurrence == "de-novo" & condition == "CTRL") %>% pull(variant)
KI36_dnV <- mtVar.bulk.df %>% filter(occurrence != "Parental" & condition == "KI36") %>% pull(variant)
KIA2_dnV <- mtVar.bulk.df %>% filter(occurrence != "Parental" & condition == "KIA2") %>% pull(variant)


mmatMut.CTRL <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.CTRL.rds")[CTRL_dnV,]
mmatMut.KI36 <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KI36.rds")[KI36_dnV,]
mmatMut.KIA2 <- readRDS("../POLG_HEK_large_data_files/output/2_mmtx.KIA2.rds")[KIA2_dnV,]
mmat.list <- list(mmatMut.CTRL, mmatMut.KI36, mmatMut.KIA2)
HEK_metadata <- read.csv("../output/1_HEK_metadata.csv",header = TRUE, row.names = 1)


## the number of mutations detected in each cell with differnt heteroplasmy threshold
# *Medium and mean number of mutation with various cutoff*
MutNum_continuous <- as.data.frame(do.call(cbind, lapply(seq(0, 1, by = 0.05), function(i){
  unlist(sapply(mmat.list, function(mmat){colSums(assays(mmat)$allele_frequency > i)}))
})))
colnames(MutNum_continuous) <- paste0("cutoff_", seq(0, 1, by = 0.05))
MutNum_continuous <- merge(MutNum_continuous, HEK_metadata[,c("condition", "atacUMAP_1", "atacUMAP_2")], by=0, all.x=TRUE) %>% tibble::column_to_rownames(var = "Row.names")

MutNum.summary <- as.data.frame(do.call(rbind, psych::describeBy(MutNum_continuous,
                                                                 group=MutNum_continuous$condition))) %>%
  tibble::rownames_to_column("condition") %>%
  filter(grepl("cutoff", condition)) %>%
  tidyr::separate_wider_delim(condition, delim = "_", names = c("condition", "cutoff")) %>%
  mutate(condition = stringr::str_replace(condition, ".cutoff", ""),
         cutoff = as.numeric(cutoff))


p1 <- MutNum.summary %>% ggplot(aes(x = cutoff, y = mean + 1, color = condition)) +
    geom_point() + geom_line() + scale_color_manual(values = color_vec) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05), add = c(0, 0))) +
    scale_y_log10(limits = c(1, 600), expand = c(0,0),
                  breaks = c(10, 100, 500),
                  labels = c(10, 100, 500)) +
    theme_classic() + labs(y = "Mean count of mutation", x = "Heteroplasmy Thresold") + Seurat::NoLegend()
ggsave(plot = p1, "../plot/41_MMB_mean_cutoff_log.pdf", width = 3, height = 3)

MutNum <- as.data.frame(do.call(cbind, lapply(c(0, 0.01, 0.03, 0.1, 0.25, 0.5, 0.75, 0.99), function(i){
  unlist(sapply(mmat.list, function(mmat){colSums(assays(mmat)$allele_frequency > i)}))
})))

colnames(MutNum) <- paste0("cut_off_", c(0, 0.01, 0.03, 0.1, 0.25, 0.5, 0.75, 0.99))

MutNum <- merge(MutNum, HEK_metadata[,c("condition", "atacUMAP_1", "atacUMAP_2", "mtDNA_depth")], by=0, all.x=TRUE) %>% 
  tibble::column_to_rownames(var = "Row.names")

# **UMAPs projected with the number of mutations at a defined heteroplasmy**
p2 <- MutNum %>% tidyr::pivot_longer(cols = -c("condition", "atacUMAP_1", "atacUMAP_2", "mtDNA_depth"), 
                                     names_to = "cut_off", values_to = "Num") %>% 
  ggplot(aes(x=atacUMAP_1, y = atacUMAP_2, col = log10(Num + 1))) +
  geom_point(size = 0.01) + # scale_colour_gradientn(colours = c("darkred", "orange", "yellow"))++ 
  scale_colour_viridis_c() +
  facet_wrap(~cut_off, ncol = 4) +
  theme_classic()
p2.nl <- p2 + theme_void() + Seurat::NoLegend() + theme(strip.background = element_blank(),strip.text.x = element_blank())
p2.leg <- ggpubr::get_legend(p2) %>% ggpubr::as_ggplot()
ggsave("../plot/42_mutNum_sample_cutoff_UMAP.pdf", p2.nl, width = 8, height = 4)
ggsave("../plot/42_mutNum_sample_cutoff_UMAP_ori.pdf", p2, width = 8, height = 4)
ggsave("../plot/42_leg.pdf", p2.leg, width = 2, height = 2)

p3 <- MutNum %>% tidyr::pivot_longer(cols = -c("condition", "atacUMAP_1", "atacUMAP_2", "mtDNA_depth"), 
                                     names_to = "cut_off", values_to = "Num") %>%  
  ggplot(aes(x=mtDNA_depth, y = Num, color = condition)) +
  geom_point(size = 0.01) + scale_color_manual(values = color_vec) + #scale_colour_gradientn(colours = c("darkred", "orange", "yellow"))+
  xlim(NA, 150) +
  facet_wrap(~cut_off, ncol = 4, scales = "free") + 
  ylab("Number of mutations") +  xlab("mtDNA depth") +
  theme_classic()
p3.nl <- p3 +Seurat::NoLegend() + theme(strip.background = element_blank(),strip.text.x = element_blank())
# ggsave(plot = p3.nl, "../plot/43_depth_MMB_cutoff.pdf", width = 8, height = 4)
# ggsave(plot = p3, "../plot/43_depth_MMB_cutoff.ori.pdf", width = 8, height = 4)


# MutNum %>% tidyr::pivot_longer(cols = -c("condition", "atacUMAP_1", "atacUMAP_2", "mtDNA_depth"), 
#                                names_to = "cut_off", values_to = "Num") %>%  
#   filter(cut_off == "cut_off_0.03" & condition == "KIA2") %>% 
#   mutate(depth_cutoff = ifelse(mtDNA_depth > 34, "> 34", "<= 34")) %>% 
#   ggplot(aes(x=mtDNA_depth, y = Num, color = depth_cutoff)) +
#   geom_point(size = 0.01) + scale_color_manual(values = c("blue", "firebrick")) + #xlim(NA, 150) +
#   # facet_wrap(~condition, ncol = 4, scales = v"free_y") + 
#   geom_smooth(method = "lm") +
#   geom_vline(xintercept = 34, linetype = "dashed") +
#   ylab("Number of mutations") +  xlab("mtDNA depth") +
#   theme_classic()

# library(chngpt)
# 
# df2 <- MutNum %>% tidyr::pivot_longer(cols = -c("condition", "atacUMAP_1", "atacUMAP_2", "mtDNA_depth"), 
#                                       names_to = "cut_off", values_to = "Num") %>%  
#   filter(cut_off == "cut_off_0.03" & condition == "KIA2") 
# 
# set.seed(123)
# fit=chngptm(formula.1=Num~1, formula.2=~mtDNA_depth, df2, type="segmented", family="gaussian", 
#             var.type="bootstrap", ci.bootstrap.size=10000)
# plot(fit)






