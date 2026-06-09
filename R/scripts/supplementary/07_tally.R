# 07_tally.R — Supplementary forest plots (all padj<0.05 features) and
# tally-plot data (Rule A, SILVA) for genus-level DA across studies.
#
# The main-figure forests that used to live here have been superseded by the
# Duvallet-style heatmaps produced by 08_heatmap.R and 09_heatmap_picrust2.R,
# which show per-study effect sizes side-by-side instead of a collapsed
# top-10 forest. This script is now the producer of:
#
#   figures/supplementary/forest_supp_<disease>_<type>.pdf  # all padj<0.05
#   figures/tally_data/tally_silva_<disease>.tsv            # (feature, up, down)
#
# Tally data feeds the tally-plot skill (Python). Rule A: a genus counts as
# "up" in a study if coef > 0 & qval < 0.25, "down" if coef < 0 & qval < 0.25.
# We keep genera appearing in >= 2 studies (>=1 for LPRD since k=1).
#
# Note on file reading: some meta-analysis TSVs have trailing NUL bytes from
# cloud-sync/Windows write truncation artifacts. read_safe() strips them.

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
})

SUPP  <- file.path(FIGURES, "supplementary")
TALLY <- file.path(FIGURES, "tally_data")
dir.create(SUPP,  recursive = TRUE, showWarnings = FALSE)
dir.create(TALLY, recursive = TRUE, showWarnings = FALSE)

# Strip embedded NULs before parsing (some files have trailing NUL padding).
read_safe <- function(path) {
  raw <- readBin(path, "raw", file.info(path)$size)
  raw <- raw[raw != as.raw(0)]
  readr::read_tsv(I(rawToChar(raw)), show_col_types = FALSE)
}

# Summary forest: one row per feature, pooled estimate + 95% CI.
summary_forest <- function(df, title) {
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() +
             labs(title = paste0(title, " (no features)")))
  }
  df <- df %>%
    mutate(
      ci_lo     = estimate - 1.96 * se,
      ci_hi     = estimate + 1.96 * se,
      direction = ifelse(estimate > 0, "UP", "DN"),
      feature   = factor(feature, levels = rev(feature))
    )
  ggplot(df, aes(x = estimate, y = feature, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c(UP = "#d4652a", DN = "#3970b0"),
                       guide = "none") +
    labs(title = title, x = "pooled log2FC (DL)", y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11, face = "bold"))
}

# ---- 1. Supplementary forest figures: all padj<0.05 features ----------------
files <- list(
  silva   = "meta_analysis__silva.tsv",
  ehomd   = "meta_analysis__ehomd.tsv",
  pathway = "meta_analysis_picrust2__pathway.tsv",
  ko      = "meta_analysis_picrust2__ko.tsv"
)

for (type in names(files)) {
  m <- read_safe(file.path(RESULTS, files[[type]]))
  sig <- m %>% dplyr::filter(!is.na(padj), padj < 0.05)
  if (nrow(sig) == 0) {
    message("No padj<0.05 features for ", type, "; skip supplementary.")
    next
  }
  for (dis in unique(sig$disease)) {
    ds <- sig %>% dplyr::filter(disease == dis) %>% dplyr::arrange(padj)
    h <- max(3, min(40, 0.25 * nrow(ds) + 2))
    p <- summary_forest(
      ds,
      sprintf("%s — %s (all padj<0.05, n=%d)", toupper(type), dis, nrow(ds))
    )
    ggsave(
      file.path(SUPP, sprintf("forest_supp_%s_%s.pdf", dis, type)),
      p, width = 8, height = h, limitsize = FALSE
    )
    message("Supp forest: ", dis, " x ", type, " (", nrow(ds), " features)")
  }
}

# ---- 2. Tally-plot data (SILVA, Rule A: qval < 0.25) ------------------------
all_silva <- read_safe(file.path(RESULTS, "maaslin2_all_studies__silva.tsv"))

for (dis in unique(all_silva$disease)) {
  dsig <- all_silva %>%
    dplyr::filter(disease == dis, !is.na(qval), qval < 0.25)

  if (nrow(dsig) == 0) {
    message("Tally (SILVA, ", dis, "): 0 sig genera under Rule A; skip.")
    next
  }

  tally <- dsig %>%
    dplyr::mutate(dir = ifelse(coef > 0, "up", "down")) %>%
    dplyr::count(feature, dir) %>%
    tidyr::pivot_wider(names_from = dir, values_from = n, values_fill = 0)

  if (!"up"   %in% names(tally)) tally$up   <- 0L
  if (!"down" %in% names(tally)) tally$down <- 0L

  min_studies <- if (dis == "LPRD") 1 else 2
  tally <- tally %>%
    dplyr::filter((up + down) >= min_studies) %>%
    dplyr::select(feature, up, down) %>%
    dplyr::arrange(desc(up - down))

  out <- file.path(TALLY, sprintf("tally_silva_%s.tsv", dis))
  readr::write_tsv(tally, out)
  message(sprintf("Tally (SILVA, %s): %d genera (>=%d studies) -> %s",
                  dis, nrow(tally), min_studies, basename(out)))
}

message("\nDone. Supplementary forests in ", SUPP)
message("Tally TSVs in ", TALLY, " — feed to the tally-plot skill.")
