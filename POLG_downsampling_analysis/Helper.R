##### make_cov_plot
make_cov_plot <- function(df, title = NULL) {
  cov_df <- as.data.frame(df, check.names = FALSE)

  if (!"pos" %in% names(df) && !is.null(rownames(df))) {
    cov_df <- cbind(pos = as.integer(rownames(df)), df)
  }

  long <- df %>%
    pivot_longer(-pos, names_to = "sample", values_to = "coverage") %>%
    mutate(
      group = str_extract(sample, "^[A-Za-z0-9]+"),   # e.g. CTRL, KIA2, KI36
      depth = str_extract(sample, "[0-9]+M")          # e.g. 10M, 20M, 50M
    ) %>% filter(depth %in% c("10M", "20M", "30M", "40M"))

  ggplot(long, aes(x = pos, y = coverage, color = group)) +
    geom_line() +
    facet_grid(depth~group, scale = "free_y") +
    scale_y_continuous(limits = c(-3, NA), expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(0, 16569, by = 5000), expand = c(0, 0)) +
    scale_color_manual(values = color_vec)+
    labs(
      title = title,
      x = "Position on mtDNA chromosome",
      y = "Roll mean coverage",
      color = "Sampling depth"
    ) +
    geom_rect(
      data = mitochondrial_data,
      aes(xmin = Starting, xmax = Ending, ymin = -3, ymax = 0, fill = color),
      inherit.aes = FALSE,
      alpha = 1
    ) +
    scale_fill_manual(values = mtGene_colors, name = "Gene") +
    guides(fill = guide_legend(ncol = 2)) +
    theme_classic()
}

### Count Mut in single cell
CountscMut <- function(x){
    df <- as.data.frame(do.call(cbind,lapply(seq(0, 1, by = 0.05), function(i){
  sapply(mmat.het.list[x], function(mmat){colSums(assays(mmat)$allele_frequency > i)})
})))

colnames(df) <- paste0(colnames(df), "_", seq(0, 1, by = 0.05) *100)
df <- df %>% tibble::rownames_to_column(var = "barcode") %>% tidyr::pivot_longer(!barcode, names_to = "label", values_to = "count") 

df <- df %>%
  extract(
    label,
    into = c("sample", "reads_in_M", "threshold"),
    regex = "^([^_]+)_(\\d+)M_(\\d+)$",
    remove = FALSE
  ) %>%
  mutate(across(c(reads_in_M, threshold), as.integer))

df   
}

### make venn diagram
# ---------- helpers ----------
get_set <- function(x, key) {
  v <- x[[key]]
  if (is.null(v)) character(0) else unique(as.character(v))
}

# Build sets of samples for a single depth (CTRL/KI36/KIA2 at, say, 10M)
get_sets_samples_at_depth <- function(var.list, depth_M, samples = c("CTRL","KI36","KIA2")) {
  sets <- setNames(
    lapply(samples, function(s) get_set(var.list, sprintf("%s_%dM", s, depth_M))),
    samples
  )
  Filter(function(z) length(z) > 0, sets) # drop empty
}

# Build sets of depths for a single sample (e.g. CTRL: 10M,20M,30M,40M)
get_sets_depths_for_sample <- function(var.list, sample, depths = c(10,20,30,40)) {
  lbls <- paste0(depths, "M")
  sets <- setNames(
    lapply(depths, function(d) get_set(var.list, sprintf("%s_%dM", sample, d))),
    lbls
  )
  Filter(function(z) length(z) > 0, sets) # drop empty
}

# ---------- plotting ----------
# Your sample colors
default_sample_colors <- c("CTRL" = "#2780FF", "KI36" = "#960096", "KIA2" = "#000090")

# Default depth colors (edit if you like)
default_depth_colors <- c("10M"="#66c2a5","20M"="#fc8d62","30M"="#8da0cb","40M"="#e78ac3")

plot_samples_at_depth <- function(var.list,
                                  depth_M,
                                  samples = c("CTRL","KI36","KIA2"),
                                  sample_colors = default_sample_colors) {
  sets <- get_sets_samples_at_depth(var.list, depth_M, samples)
  if (length(sets) < 2L) warning(sprintf("Not enough groups at %dM to draw a Venn.", depth_M))

  fill_cols <- unname(sample_colors[names(sets)])
  fill_cols[is.na(fill_cols)] <- "#BBBBBB"

  ggvenn(
    sets,
    fill_color     = fill_cols,
    fill_alpha     = 0.45,
    stroke_color   = "grey30",
    stroke_size    = 0.6,
    show_elements  = FALSE,
    show_percentage= FALSE,
    set_name_size  = 5,
    text_size      = 4
  ) +
    ggtitle(sprintf("Samples @ %dM", depth_M)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

plot_depths_within_sample <- function(var.list,
                                      sample,
                                      depths = c(10,20,30,40),
                                      depth_colors = default_depth_colors) {
  sets <- get_sets_depths_for_sample(var.list, sample, depths)
  if (length(sets) < 2L) warning(sprintf("Not enough depths for %s to draw a Venn.", sample))

  fill_cols <- unname(depth_colors[names(sets)])
  fill_cols[is.na(fill_cols)] <- "#BBBBBB"

  ggvenn(
    sets,
    fill_color     = fill_cols,
    fill_alpha     = 0.45,
    stroke_color   = "grey30",
    stroke_size    = 0.6,
    show_elements  = FALSE,
    show_percentage= FALSE,
    set_name_size  = 5,
    text_size      = 4
  ) +
    ggtitle(sprintf("%s: overlaps across depths (%s)",
                    sample, paste(names(sets), collapse = ", "))) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# ---------- one-click exporters ----------
make_venn_pdfs <- function(var.list,
                           depths  = c(10,20,30,40),
                           samples = c("CTRL","KI36","KIA2"),
                           sample_colors = default_sample_colors,
                           depth_colors  = default_depth_colors,
                           out_prefix = "venn") {

  # 1) Samples-at-depth: one page per depth
  pdf(sprintf("%s_samples-at-depth_%s.pdf",
              out_prefix, paste0(depths, "M", collapse = "-")), width = 7, height = 6)
  for (d in depths) print(plot_samples_at_depth(var.list, d, samples, sample_colors))
  dev.off()

  # 2) Depths-within-sample: one page per sample
  pdf(sprintf("%s_depths-within-sample_%s.pdf",
              out_prefix, paste0(depths, "M", collapse = "-")), width = 7, height = 6)
  for (s in samples) print(plot_depths_within_sample(var.list, s, depths, depth_colors))
  dev.off()
}

### get only non-parental variants
# --- helpers ---
get_sample <- function(k) sub("_.*$", "", k)
get_depth  <- function(k) as.integer(sub("^.*_([0-9]+)M$", "\\1", k))
get_key    <- function(s, d) sprintf("%s_%dM", s, d)

get_variants <- function(lst, s, d) {
  v <- lst[[get_key(s, d)]]
  v <- if (is.null(v)) character(0) else unique(as.character(v))
  v[!is.na(v)]
}

# union of a list of character vectors
u_union <- function(xlist) {
  if (length(xlist) == 0) character(0) else unique(unlist(xlist, use.names = FALSE))
}

# fallback control set: CTRL_d if exists; else CTRL_40M; else union of CTRL_10/20/30/40; else empty
ctrl_ref <- function(lst, d, fallback_depth = 40L) {
  v <- get_variants(lst, "CTRL", d)
  if (length(v) > 0) return(v)
  v <- get_variants(lst, "CTRL", fallback_depth)
  if (length(v) > 0) return(v)
  v <- u_union(lapply(c(10L,20L,30L,40L), function(dd) get_variants(lst, "CTRL", dd)))
  v
}

# --- main filter ---
filter_variants_relative_to_CTRL <- function(var.list, fallback_ctrl_depth = 40L,
                                             samples = c("CTRL","KI36","KIA2")) {
  keys <- names(var.list)
  out  <- vector("list", length(keys)); names(out) <- keys
  summ <- data.frame(sample=character(), depth=integer(),
                     n_before=integer(), n_after=integer(),
                     stringsAsFactors = FALSE)

  # pre-compute per-depth unions for non-CTRL samples (used for CTRL unique-at-depth)
  all_depths <- sort(unique(get_depth(keys)))
  other_map  <- lapply(all_depths, function(d) {
    u_union(lapply(setdiff(samples, "CTRL"), function(s) get_variants(var.list, s, d)))
  })
  names(other_map) <- as.character(all_depths)

  for (k in keys) {
    s <- get_sample(k)
    d <- get_depth(k)
    v <- get_variants(var.list, s, d)

    if (s == "CTRL") {
      # CTRL_d minus union(KI36_d, KIA2_d)   (only at that depth if present)
      others_d <- other_map[[as.character(d)]]
      ref <- others_d
      filtered <- if (length(ref) > 0) setdiff(v, ref) else v
    } else {
      # Non-CTRL: sample_d minus (CTRL_d if exists, else CTRL_40M)
      ref <- ctrl_ref(var.list, d, fallback_ctrl_depth)
      filtered <- setdiff(v, ref)
    }

    out[[k]] <- filtered
    summ <- rbind(summ, data.frame(sample=s, depth=d,
                                   n_before=length(v), n_after=length(filtered)))
  }

  attr(out, "summary") <- summ[order(summ$sample, summ$depth), ]
  out
}
