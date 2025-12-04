mgatk_filter <- function(mmat, n_cells = 1, vmr.thres = 0.005){
  var <- as.data.frame(rowData(mmat))
  var <- var %>% subset(n_cells_conf_detected >= n_cells &
                          strand_correlation > 0.65 &
                          vmr >= vmr.thres &
                          variant %ni% c("301A>C", "302A>C", "310T>C", "316G>C")) %>% pull(variant)
  mmat <- mmat[var,]
  return(mmat)
}

# MITOMAP disease variants
mitomap.disease <- read.table(
  "MITOMAP_disease_2022-05-25.txt",
  header = TRUE, sep = "\t", na.strings = c("", "NA"),
  fill = TRUE, comment.char = "", quote = ""
)
# MLC score
MLC <- read.table("../../global/MLC_score.tsv", header = TRUE)

# mitoVAR annotation
mito.annotation <- read.table("mitovar_annotation.txt", header = TRUE) %>%
  mutate(
    variant = gsub("MT_(\\d+)_([A-Z])/([A-Z])", "\\1\\2>\\3", Uploaded_variation),
    Position = as.integer(sub("([0-9]+).*", "\\1", variant))
  ) %>%
  left_join(MLC, by = "Position") %>%
  select(variant, Consequence, SYMBOL, BIOTYPE, Protein_position, Amino_acids, MLC_pos_score, SIFT, PolyPhen) %>%
  mutate(
    effect = case_when(
      Consequence %in% c("stop_retained_variant", "synonymous_variant") ~ "synonymous",
      Consequence %in% c("missense_variant") ~ "missense",
      Consequence %in% c("intergenic_variant", "non_coding_transcript_exon_variant") ~ "non-coding",
      TRUE ~ "truncating"),
    SIFT_score = ifelse(
      SIFT == "-",
      NA,
      as.numeric(sub(".*\\(([^)]+)\\)", "\\1", SIFT))),
    PolyPhen_score = ifelse(
      PolyPhen == "-",
      NA,
      as.numeric(sub(".*\\(([^)]+)\\)", "\\1", PolyPhen))),
  )

# Define confirmed pathogenic variants
disease.var <- mitomap.disease %>%
  mutate(variant = paste0(pos, ref, ">", alt)) %>%
  filter(status %in% c("Cfrm", "Cfrm [P*]", " Cfrm [LP*]")) %>%
  pull(variant)

# Tag pathogenic
mito.annotation <- mito.annotation %>%
  mutate(effect = ifelse(variant %in% disease.var, "pathogenic", effect))
