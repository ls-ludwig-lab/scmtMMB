setwd("~/Ludwig_lab/scmtMMB/PBMC/code/")
### Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(data.table)
library(ggpubr)
library(purrr)
'%ni%' <- Negate('%in%')
## Calculate MMB
color_vec <- c("M80" = "darkorchid4","M60" = "firebrick4", "M35" = "violetred3", "M29" = "firebrick1", "H47" = "blue3", "H05" = "dodgerblue3")

# load metadata
scmtMMB <- read.csv("../PBMC_large_data_files/output/2_scmtMMB.csv", row.names = "X")
scmtMMB_ex <- read.csv("../PBMC_large_data_files/output/2_scmtMMB_exclude_3243.csv", row.names = "X")

## mtscMMB cor with depth: correlation plotting
theme_cor <- theme(aspect.ratio=1/1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
corplot.MPM <- scmtMMB %>% dplyr::filter(symbol == "Total") %>%
  ggplot(aes(x = log10(mtDNA_depth), y = mutation_per_MB, col = sample)) +
  geom_point(size = 0.01) + facet_wrap(~sample, ncol = 6) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("Total scMPM") + theme(legend.position="none") +
  geom_smooth(method = "lm") + theme_cor

ggsave(plot = corplot.MMB, "../plot/21_cor_depth_MPM.pdf", width=9, height=2)


corplot.MSS <- scmtMMB_ex %>% dplyr::filter(symbol == "Total") %>%
  ggplot(aes(x = log10(mtDNA_depth), y = MSS_weighted, col = sample)) +
  geom_point(size = 0.01) + facet_wrap(~sample, ncol = 6) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("Total scwMSS") + theme(legend.position="none") +
  geom_smooth(method = "lm") + theme_cor
ggsave(plot = corplot.MSS, "../plot/21_cor_depth_MSS.pdf", width=9, height=2)

corplot.MMB <- scmtMMB_ex %>% dplyr::filter(symbol == "Total") %>%
  ggplot(aes(x = mutation_per_MB, y = MSS_weighted, col = sample)) +
  geom_point(size = 0.01) + facet_wrap(~sample, ncol = 6) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("Total scwMSS") + theme(legend.position="none") + xlab("Total scMPM") +
  geom_smooth(method = "lm") + theme_cor
ggsave(plot = corplot.MMB, "../plot/21_cor_MPM_MSS.pdf", width=9, height=2)

## mtscMMB cor with depth: linear regression


all_regress <- scmtMMB %>% 
  filter(symbol == "Total") %>% 
  group_by(sample) %>% 
  nest() %>% 
  mutate(
    #mod1 = map(data, ~ lm(mutation_per_MB ~ log10(mtDNA_depth), data = .x)),
    rp.MPM = map_dbl(data, ~ cor(.x$mutation_per_MB, log10(.x$mtDNA_depth), method = "pearson", use = "complete.obs")),
    rs.MPM = map_dbl(data, ~ cor(.x$mutation_per_MB, log10(.x$mtDNA_depth), method = "spearman", use = "complete.obs")),
    #mod2 = map(data, ~ lm(MSS_weighted ~ log10(mtDNA_depth), data = .x)),
    rp.MSS = map_dbl(data, ~ cor(.x$MSS_weighted, log10(.x$mtDNA_depth), method = "pearson", use = "complete.obs")),
    rs.MSS = map_dbl(data, ~ cor(.x$MSS_weighted, log10(.x$mtDNA_depth), method = "spearman", use = "complete.obs")),
    #mod3 = map(data, ~ lm(MSS_weighted ~ mutation_per_MB, data = .x)),
    rp.MMB = map_dbl(data, ~ cor(.x$MSS_weighted, .x$mutation_per_MB, method = "pearson", use = "complete.obs")),
    rs.MMB = map_dbl(data, ~ cor(.x$MSS_weighted, .x$mutation_per_MB, method = "spearman", use = "complete.obs"))
  ) %>% 
  ungroup()

# A tibble: 6 × 8
# sample data                    rp.MPM  rs.MPM  rp.MSS  rs.MSS rp.MMB rs.MMB
# <chr>  <list>                   <dbl>   <dbl>   <dbl>   <dbl>  <dbl>  <dbl>
# 1 M80    <tibble [5,060 × 16]>  0.0848   0.0704  0.0674  0.0412  0.891  0.906
# 2 M60    <tibble [6,223 × 16]>  0.00963 -0.0150 -0.0745 -0.0998  0.815  0.849
# 3 M35    <tibble [6,412 × 16]>  0.102    0.120   0.0980  0.111   0.910  0.925
# 4 M29    <tibble [5,743 × 16]>  0.183    0.218   0.197   0.230   0.925  0.941
# 5 H47    <tibble [6,204 × 16]>  0.0942   0.146   0.100   0.144   0.850  0.891
# 6 H05    <tibble [15,662 × 16]> 0.0448   0.0611  0.0401  0.0733  0.845  0.882

#use broom the extract the slope and rsq per group
# glance <-all_regress %>% 
#   mutate(tidy1 = map(mod1, broom::tidy), glance1 = map(mod1, broom::glance), augment1 = map(mod1, broom::augment),
#          rsq.MPM = glance1 %>% map_dbl('r.squared'), pval.MPM = glance1 %>% map_dbl('p.value'), slope.MPM = tidy1 %>% map_dbl(function(x) x$estimate[2]),
#          tidy2 = map(mod2, broom::tidy), glance2 = map(mod2, broom::glance), augment2 = map(mod2, broom::augment),
#          rsq.MSS = glance2 %>% map_dbl('r.squared'), pval.MSS = glance2 %>% map_dbl('p.value'), slope.MSS = tidy2 %>% map_dbl(function(x) x$estimate[2]),                    
#          tidy3 = map(mod3, broom::tidy), glance3 = map(mod3, broom::glance), augment3 = map(mod3, broom::augment),
#          rsq.MMB = glance3 %>% map_dbl('r.squared'), pval.MMB = glance3 %>% map_dbl('p.value'), slope.MMB = tidy3 %>% map_dbl(function(x) x$estimate[2])) 
# glance

scmtMMB_ex %>% mutate(age = as.numeric(gsub("H|M", "", sample))) %>% filter(symbol == "Total") %>% 
  ggplot(aes(x= age, y = mutation_per_MB)) +
  geom_violin(aes(fill = sample), scale = "width") + scale_fill_manual(values = color_vec) + 
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), shape = 23, size = 2, fill = "white") +
  Seurat::NoLegend() + geom_smooth(method = "lm") + theme_bw() + 
  theme(aspect.ratio=1/1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

model = lm(mutation_per_MB ~age, data = scmtMMB_ex %>% mutate(age = as.numeric(gsub("H|M", "", sample))) %>% filter(symbol == "Total"))
mutate(tidy = map(mod, broom::tidy), glance = map(mod, broom::glance), augment1 = map(mod, broom::augment),
       rsq.MPM = glance1 %>% map_dbl('r.squared'), pval.MPM = glance1 %>% map_dbl('p.value'), slope.MPM = tidy1 %>% map_dbl(function(x) x$estimate[2]))

### age-dependent increase of MMB
scmtMMB %>% mutate(age = as.numeric(gsub("H|M", "", sample))) %>% filter(symbol == "Total") %>% 
  ggplot(aes(x= age, y = mutation_per_MB)) +
  geom_violin(aes(fill = sample), scale = "width") + scale_fill_manual(values = color_vec) + 
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), shape = 23, size = 1, fill = "white") +
  Seurat::NoLegend() + geom_smooth(method = "lm") + theme_bw() + 
  theme(aspect.ratio=1/1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + Seurat::NoLegend()

scmtMMB_ex %>% mutate(age = as.numeric(gsub("H|M", "", sample))) %>% filter(symbol == "Total") %>% 
  ggplot(aes(x= age, y = mutation_per_MB)) +
  geom_violin(aes(fill = sample), scale = "width") + scale_fill_manual(values = color_vec) + 
  geom_point(data = . %>% dplyr::filter(symbol == "Total" & barcode == "pseudobulk"), shape = 23, size = 1, fill = "white") +
  Seurat::NoLegend() + geom_smooth(method = "lm") + theme_bw() + 
  theme(aspect.ratio=1/1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Total scMPM") + xlab("Age (Years)") + Seurat::NoLegend()

model = lm(mutation_per_MB ~age, 
           data = scmtMMB_ex %>% mutate(age = as.numeric(gsub("H|M", "", sample))) %>% filter(symbol == "Total"))
# Coefficients:
# (Intercept)          age  
#     29.7594       0.7425  
broom::glance(model)
# r.squared adj.r.squared sigma statistic p.value    df   logLik     AIC     BIC  deviance df.residual  nobs
#     <dbl>         <dbl> <dbl>     <dbl>   <dbl> <dbl>    <dbl>   <dbl>   <dbl>     <dbl>       <int> <int>
#     0.253         0.253  32.4    15373.       0     1 -221869. 443743. 443770. 47580646.       45302 45304

scmtMMB %>% filter(grepl(x = sample, "M")) %>% mutate(age = as.numeric(gsub("H|M", "", sample))) %>% 
  filter(symbol == "Total") %>% ggplot(aes(x = X3243A_G, colour = sample)) + 
  stat_ecdf() + scale_color_manual(values = color_vec) + theme_bw() +
  theme(aspect.ratio=1/1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

scmtMMB %>% filter(grepl(x = sample, "M")) %>% filter(symbol == "Total") %>% 
  ggplot(aes(x= sample, y = X3243A_G)) + geom_jitter(color = "grey", size =0.01) +
  geom_violin(aes(fill = sample), scale = "width") + 
  scale_fill_manual(values = color_vec) + 
  Seurat::NoLegend() + geom_smooth(method = "lm") + theme_bw() + 
  theme(aspect.ratio=1/1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + Seurat::NoLegend()

scmtMMB %>% filter(grepl(x = sample, "M")) %>% filter(symbol == "tRNA") %>% 
  ggplot(aes(x= X3243A_G, y = mutation_per_MB, color = sample)) + geom_point() + geom_contour(aes(z = mtDNA_depth))+
  scale_color_manual(values = color_vec) + 
  Seurat::NoLegend() + geom_smooth(method = "lm", color = "black") + theme_bw() + facet_wrap(~sample) +
  theme(aspect.ratio=1/1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + Seurat::NoLegend()

## mtscMMB cor with m.3242A>G
theme_cor <- theme(aspect.ratio=1/1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
corplot.MPM_3243 <- scmtMMB %>% dplyr::filter(symbol == "Total") %>% filter(grepl(x = sample, "M")) %>% 
  ggplot(aes(x = X3243A_G * 100, y = mutation_per_MB, col = sample)) +
  geom_point(size = 0.01) + facet_wrap(~sample, ncol = 4) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("Total scmtMPM")+ xlab("m.3243A>G heteroplasmy (%)") + theme(legend.position="none") +
  geom_smooth(method = "lm", color = "black") + theme_cor
ggsave(plot = corplot.MPM_3243, "../plot/21_cor_3243_scmtMPM.pdf", width=4, height=2)

corplot.MPM_3243_density <- scmtMMB %>% dplyr::filter(symbol == "Total") %>% filter(sample %in% c("M29","M35")) %>% 
  ggplot(aes(x = X3243A_G * 100, y = mutation_per_MB, col = sample)) + 
  # geom_point(color = "grey75", size = 0.01) +
  geom_density2d(bins = 50) + scale_x_continuous(limits = c(-18, 110), breaks = seq(0,100, 25))+
  scale_y_continuous(limits = c(-18, 160), breaks = seq(0,150, 50))+
  facet_wrap(~sample, ncol = 6) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("Total scmtMPM")+ xlab("m.3243A>G heteroplasmy (%)") + theme(legend.position="none") +
  theme_cor

ggsave(plot = corplot.MPM_3243_density, "../plot/21_cor_3243_MPM_d.pdf", width=4, height=2)


corplot.MPM_3243_density <- scmtMMB_ex %>% dplyr::filter(symbol == "Total") %>% filter(sample %in% c("M29","M35")) %>% 
  ggplot(aes(x = X3243A_G * 100, y = mutation_per_MB, col = sample)) + 
  # geom_point(color = "grey75", size = 0.01) +
  geom_density2d(bins = 50) + scale_x_continuous(limits = c(-18, 110), breaks = seq(0,100, 25))+
  scale_y_continuous(limits = c(-18, 160), breaks = seq(0,150, 50))+
  facet_wrap(~sample, ncol = 6) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("Total scmtMPM")+ xlab("m.3243A>G heteroplasmy (%)") + theme(legend.position="none") +
  theme_cor

ggsave(plot = corplot.MPM_3243_density, "../plot/21_cor_ex3243_MPM_d.pdf", width=4, height=2)

#scwMSS vs 3243 VAF
corplot.MSS_3243_density <- scmtMMB %>% dplyr::filter(symbol == "Total") %>% filter(sample %in% c("M29","M35")) %>% 
  ggplot(aes(x = X3243A_G * 100, y = MSS_weighted, col = sample)) + 
  #geom_point(color = "grey75", size = 0.01) +
  geom_density2d(bins = 100) + scale_x_continuous(limits = c(-18, 110), breaks = seq(0,100, 25))+
  facet_wrap(~sample, ncol = 6) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("Total scwMSS")+ xlab("m.3243A>G heteroplasmy (%)") + theme(legend.position="none") +
  theme_cor

ggsave(plot = corplot.MSS_3243_density, "../plot/21_cor_3243_MSS_d.pdf", width=4, height=2)


corplot.depth_MPM <- scmtMMB %>% dplyr::filter(symbol == "Total") %>% 
  filter(sample %in% c("M29","M35")) %>% filter(X3243A_G == 0) %>% 
  ggplot(aes(x = mtDNA_depth, y = mutation_per_MB, col = sample)) + 
  geom_point(color = "black", size = 0.01) +
  # geom_density2d(bins = 100) + scale_x_continuous(limits = c(-18, 110), breaks = seq(0,100, 25))+
  facet_wrap(~sample, ncol = 6) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("scmtMPM")+ xlab("mtDNA depth") + theme(legend.position="none") +
  theme_cor
ggsave(plot = corplot.depth_MPM, "../plot/21_cor_MPM_depth_3243_0.pdf", width=4, height=2)


scmtMMB %>% dplyr::filter(symbol == "Total") %>% 
  filter(sample %in% c("M29","M35")) %>% filter(X3243A_G == 0) %>% 
  ggplot(aes(x = mtDNA_depth, y = mutation_per_MB, col = sample)) + 
  geom_point(color = "black", size = 0.01) +
  # geom_density2d(bins = 100) + scale_x_continuous(limits = c(-18, 110), breaks = seq(0,100, 25))+
  facet_wrap(~sample, ncol = 6) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("scmtMPM")+ xlab("mtDNA depth") + theme(legend.position="none") +
  theme_cor

VAF_depth <- scmtMMB %>% dplyr::filter(symbol == "Total") %>% 
  filter(sample %in% c("M29","M35")) %>% # filter(X3243A_G == 0) %>% 
  ggplot(aes(x = X3243A_G, y = mtDNA_depth, col = sample)) + 
  geom_point(size = 0.01) +
  # geom_density2d(bins = 100) + scale_x_continuous(limits = c(-18, 110), breaks = seq(0,100, 25))+
  facet_wrap(~sample, ncol = 2) +  scale_color_manual(values = color_vec) +
  theme_bw() + ylab("mtDNA depth")+ xlab("m.3243A>G heteroplasmy (%)") + theme(legend.position="none") +
  theme_cor
ggsave(plot = VAF_depth, "../plot/21_depth_3243.pdf", width=4, height=2)

# tRNA score vs heteroplasmy (sanity check)
corplot.tRNA_MPM_3243 <- scmtMMB %>% dplyr::filter(symbol == "tRNA" & sample %in% c("M29","M35", "M60","M80")) %>%
  ggplot(aes(x = X3243A_G * 100, y = mutation_per_MB, col = sample)) +
  geom_point(size = 0.01) + facet_wrap(~sample, ncol = 4) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("tRNA scmtMPM")+ xlab("m.3243A>G heteroplasmy (%)") + theme(legend.position="none") +
  geom_smooth(method = "lm", color = "black") + theme_cor
ggsave(plot = corplot.tRNA_MPM_3243, "../plot/21_cor_3243_tRNA_scmtMPM.pdf", width=4, height=2)



scmtMMB %>% dplyr::filter(symbol == "tRNA" & sample %in% c("M29","M35")) %>%
  ggplot(aes(x = X3243A_G * 100, y = MSS_weighted, col = sample)) +
  geom_point(size = 0.01) + facet_wrap(~sample, ncol = 6) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("tRNA scwMSS")+ xlab("m.3243A>G heteroplasmy (%)") + theme(legend.position="none") +
  geom_smooth(method = "lm", color = "black") + theme_cor

tRNA_regress_3243 <-  scmtMMB  %>% dplyr::filter(symbol == "tRNA" & sample %in% c("M29","M35", "M60","M80")) %>%
  mutate(X3243A_G = X3243A_G*100) %>% group_by(sample) %>%
  do(mod1 = lm(mutation_per_MB ~ X3243A_G, data = .),
     mod2 = lm(MSS_weighted ~ X3243A_G, data = .)) %>% ungroup() 

scmtMMB %>% dplyr::filter(symbol %in% c("Total","tRNA") & sample %in% c("M29","M35", "M60","M80")) %>%
  mutate(X3243A_G = X3243A_G*100) %>% group_by(symbol,sample) %>% 
  nest() %>% 
  mutate(
    #mod1 = map(data, ~ lm(mutation_per_MB ~ log10(mtDNA_depth), data = .x)),
    rp.MPM = map_dbl(data, ~ cor(.x$mutation_per_MB, .x$X3243A_G, method = "pearson", use = "complete.obs")),
    rs.MPM = map_dbl(data, ~ cor(.x$mutation_per_MB, .x$X3243A_G, method = "spearman", use = "complete.obs")),
    #mod2 = map(data, ~ lm(MSS_weighted ~ log10(mtDNA_depth), data = .x)),
    rp.MSS = map_dbl(data, ~ cor(.x$MSS_weighted, .x$X3243A_G, method = "pearson", use = "complete.obs")),
    rs.MSS = map_dbl(data, ~ cor(.x$MSS_weighted, .x$X3243A_G, method = "spearman", use = "complete.obs"))
    ) %>% 
  ungroup()

# # A tibble: 8 × 7
#   sample symbol data                  rp.MPM rs.MPM rp.MSS rs.MSS
#   <chr>  <chr>  <list>                 <dbl>  <dbl>  <dbl>  <dbl>
# 1 M80    Total  <tibble [5,060 × 15]>  0.260  0.178  0.305  0.219
# 2 M80    tRNA   <tibble [5,060 × 15]>  0.806  0.538  0.877  0.570
# 3 M60    Total  <tibble [6,223 × 15]>  0.218  0.160  0.144  0.102
# 4 M60    tRNA   <tibble [6,223 × 15]>  0.787  0.403  0.855  0.418
# 5 M35    Total  <tibble [6,412 × 15]>  0.648  0.638  0.734  0.709
# 6 M35    tRNA   <tibble [6,412 × 15]>  0.934  0.864  0.976  0.883
# 7 M29    Total  <tibble [5,743 × 15]>  0.643  0.608  0.729  0.668
# 8 M29    tRNA   <tibble [5,743 × 15]>  0.944  0.814  0.974  0.826


tRNA_glance_3243 <-tRNA_regress_3243 %>% 
  mutate(tidy1 = map(mod1, broom::tidy), glance1 = map(mod1, broom::glance), augment1 = map(mod1, broom::augment),
         rsq.MPM = glance1 %>% map_dbl('r.squared'), pval.MPM = glance1 %>% map_dbl('p.value'), slope.MPM = tidy1 %>% map_dbl(function(x) x$estimate[2]),
         tidy2 = map(mod2, broom::tidy), glance2 = map(mod2, broom::glance), augment2 = map(mod2, broom::augment),
         rsq.MSS = glance2 %>% map_dbl('r.squared'), pval.MSS = glance2 %>% map_dbl('p.value'), slope.MSS = tidy2 %>% map_dbl(function(x) x$estimate[2])) 
tRNA_glance_3243
# sample mod1   mod2   tidy1    glance1  augment1 rsq.MPM pval.MPM slope.MPM tidy2    glance2  augment2 rsq.MSS pval.MSS
# <chr>  <list> <list> <list>   <list>   <list>     <dbl>    <dbl>     <dbl> <list>   <list>   <list>     <dbl>    <dbl>
# 1 M29    <lm>   <lm>   <tibble> <tibble> <tibble>   0.892        0      8.63 <tibble> <tibble> <tibble>   0.949        0
# 2 M35    <lm>   <lm>   <tibble> <tibble> <tibble>   0.872        0      8.85 <tibble> <tibble> <tibble>   0.953        0
# 3 M60    <lm>   <lm>   <tibble> <tibble> <tibble>   0.620        0      8.91 <tibble> <tibble> <tibble>   0.731        0
# 4 M80    <lm>   <lm>   <tibble> <tibble> <tibble>   0.649        0      8.40 <tibble> <tibble> <tibble>   0.768        0


corplot.MPMex_3243 <- scmtMMB_ex %>% dplyr::filter(symbol == "Total") %>%
  ggplot(aes(x = X3243A_G * 100, y = mutation_per_MB, col = sample)) +
  geom_point(size = 0.01) + facet_wrap(~sample, ncol = 6) + scale_color_manual(values = color_vec) +
  theme_bw() + ylab("Total scMPM")+ xlab("m.3243A>G heteroplasmy (%)") + theme(legend.position="none") +
  geom_smooth(method = "lm", color = "black") + theme_cor
ggsave(plot = corplot.MPM_3243, "../plot/21_cor_3243_MPMex.pdf", width=10, height=2)

all_regress_3243 <-  scmtMMB %>% dplyr::filter(symbol == "Total") %>% 
  mutate(X3243A_G = X3243A_G*100) %>% group_by(sample) %>%
  do(mod1 = lm(mutation_per_MB ~ X3243A_G, data = .),
     mod2 = lm(MSS_weighted ~ X3243A_G, data = .)) %>% ungroup()

#use broom the extract the slope and rsq per group
glance_3243 <-all_regress_3243 %>% 
  mutate(tidy1 = map(mod1, broom::tidy), glance1 = map(mod1, broom::glance), augment1 = map(mod1, broom::augment),
         rsq.MPM = glance1 %>% map_dbl('r.squared'), pval.MPM = glance1 %>% map_dbl('p.value'), slope.MPM = tidy1 %>% map_dbl(function(x) x$estimate[2]),
         tidy2 = map(mod2, broom::tidy), glance2 = map(mod2, broom::glance), augment2 = map(mod2, broom::augment),
         rsq.MSS = glance2 %>% map_dbl('r.squared'), pval.MSS = glance2 %>% map_dbl('p.value'), slope.MSS = tidy2 %>% map_dbl(function(x) x$estimate[2])) 
glance_3243

all_regress_3243ex <-  scmtMMB_ex %>% dplyr::filter(symbol == "Total") %>%  
  mutate(X3243A_G = X3243A_G*100) %>% group_by(sample) %>%
  do(mod1 = lm(mutation_per_MB ~ X3243A_G, data = .),
     mod2 = lm(MSS_weighted ~ X3243A_G, data = .)) %>% ungroup()

#use broom the extract the slope and rsq per group
glance_3243ex <-all_regress_3243ex %>% 
  mutate(tidy1 = map(mod1, broom::tidy), glance1 = map(mod1, broom::glance), augment1 = map(mod1, broom::augment),
         rsq.MPM = glance1 %>% map_dbl('r.squared'), pval.MPM = glance1 %>% map_dbl('p.value'), slope.MPM = tidy1 %>% map_dbl(function(x) x$estimate[2]),
         tidy2 = map(mod2, broom::tidy), glance2 = map(mod2, broom::glance), augment2 = map(mod2, broom::augment),
         rsq.MSS = glance2 %>% map_dbl('r.squared'), pval.MSS = glance2 %>% map_dbl('p.value'), slope.MSS = tidy2 %>% map_dbl(function(x) x$estimate[2])) 
glance_3243ex



scmtMMB %>% dplyr::filter(symbol == "Total") %>% 
  filter(sample %in% c("M29","M35")) %>% filter(X3243A_G == 0) %>% 
  ggplot(aes(x = predicted.celltype.l2, y = mtDNA_depth)) + 
  geom_violin()+
  facet_wrap(~sample, ncol = 2) +  scale_color_manual(values = color_vec) +
  theme_bw() + ylab("mtDNA depth")+ xlab("m.3243A>G heteroplasmy (%)") + theme(legend.position="none") +
  theme_cor


scmtMMB %>% filter(barcode != "pseudobulk" & symbol == "Total") %>%
  mutate(
    celltype = case_when(
      predicted.celltype.l2 %in% c("CD8 TEM", "CD8 TCM") ~ "CD8 Eff/Mem",
      predicted.celltype.l2 %in% c("CD4 TEM", "CD4 TCM", "CD4 CTL") ~ "CD4 Eff/Mem",
      predicted.celltype.l2 %in% c("B naive", "B intermediate", "B memory", "Plasmablast") ~ "B cell",
      predicted.celltype.l2 %in% c("CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC") ~ "Monocyte/DC",
      grepl("NK|ILC", predicted.celltype.l2) ~ "NK/ILC",
      predicted.celltype.l2 %in% c("Doublet", "Eryth", "Platelet", "ASDC", NA) ~ "discard",
      TRUE ~ predicted.celltype.l2
    )
  ) %>% filter(celltype != "discard") %>% 
  dplyr::filter(symbol == "Total") %>% 
  filter(sample %in% c("M29","M35")) %>% filter(X3243A_G == 0) %>% 
  ggplot(aes(x = celltype, y = mtDNA_depth)) + 
  geom_violin()+
  geom_jitter(size = 0.1, color = "grey")+
  facet_wrap(~sample, ncol = 2) + theme_classic()

