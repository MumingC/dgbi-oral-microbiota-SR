# 08_heatmap.R — Duvallet-style effect-size heatmap for genus-level taxa.
# Main figure candidate (SILVA) + supplementary (eHOMD).
#
# Two configurations are emitted per database and per row-order:
#
#   all       — every genus present in meta_analysis__{db}.tsv (>=2 studies in
#               at least one disease). Meta-cell indicator: bold black border
#               when meta padj < 0.05.
#   v2        — concordance-filtered: keep genera where, within at least one
#               panel, >=2 disease-matched studies share coef sign AND >=1 of
#               those agreeing studies has qval < 0.25. Meta-cell indicator:
#               ** only at meta padj < 0.05 (no * tier). Per-study cells still
#               use the * / ** scheme at qval < 0.25 / < 0.05 because per-study
#               testing is exploratory; the Meta column is the primary meta-
#               analytic inference and uses a stricter threshold.
#
# v2 is the directional-consistency main figure. Earlier versions used a
# directional triangle on every concordant Meta cell, but because the current
# re-analysis produced 0 genera at meta padj < 0.05, triangles appeared on
# every IBS-D Meta row and visually impersonated significance markers. We now
# use a stricter ** padj < 0.05 mark on Meta cells (no * tier), so cells with
# no asterisk are concordant by inclusion criterion but not meta-significant.
#
#   Upper-GI panel: Ziganshina_2020, Qian_2023, Hao_2022, Kawar_2021,
#                   Zheng_2024, GERD-Meta
#   IBS-D    panel: Li_2025, Tang_2023, IBS-D-Meta
#
# LPRD (Zheng_2024) has k=1, so it contributes a study column but no meta
# column. It also does NOT contribute to GERD concordance.
#
# Cell annotations:
#   Per-study cells: "*" qval < STAR_THRESH_1, "**" qval < STAR_THRESH_2
#   Meta cells:      border (all) OR triangle (v2), see above.
#
# Outputs (figures/heatmap/):
#   heatmap_silva_by_padj.{pdf,png,tiff}     # supp — all features, padj border
#   heatmap_silva_by_clust.{pdf,png,tiff}    # supp — all features, padj border
#   heatmap_silva_v2_by_padj.{pdf,png,tiff}  # main — concordance features, triangle
#   heatmap_silva_v2_by_clust.{pdf,png,tiff} # main — concordance features, triangle
#   heatmap_ehomd_*                          # same pattern for eHOMD
# PNG  = 300 DPI Cairo (review / draft inclusion)
# TIFF = 600 DPI Cairo, LZW-compressed (journal submission)
#
# Reference: Duvallet et al. 2017, Nat Commun 8:1784 (Fig 2).

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))
source(here::here("lib", "heatmap_helpers.R"))

install_heatmap_deps()

HEATDIR <- file.path(FIGURES, "heatmap")
dir.create(HEATDIR, recursive = TRUE, showWarnings = FALSE)

# All-meta row filter: keep every genus that made it into the meta table.
# Accepts per_study_df for API compatibility with row_filter_concordance().
row_filter_all_meta <- function(meta_df, per_study_df = NULL) {
  meta_df %>% dplyr::distinct(feature) %>% dplyr::pull(feature)
}

# Row labels: strip "g__" prefix for readability.
strip_genus_prefix <- function(x) sub("^g__", "", x)

thresholds <- list(
  star1       = 0.25,   # per-study *  threshold
  star2       = 0.05,   # per-study ** threshold
  meta_star1  = NA,     # NA => no * tier on Meta cells (suppresses padj < 0.25)
  meta_star2  = 0.05,   # Meta ** at padj < 0.05 only
  meta_sig    = 0.05,   # border indicator threshold (used by config "border")
  color_limit = 3
)

configs <- list(
  list(tag = "",    filter_fn = row_filter_all_meta,    indicator = "border",
       desc = "all meta features, padj border"),
  list(tag = "_v2", filter_fn = row_filter_concordance, indicator = "star",
       desc = "concordance filter, meta padj ** at < 0.05 only (no * tier)")
)

for (db in c("silva", "ehomd")) {
  for (cfg in configs) {
    message("\n=== Heatmap: ", db, " [", cfg$desc, "] ===")
    dat <- build_heatmap_data(
      per_study_path = file.path(RESULTS, paste0("maaslin2_all_studies__", db, ".tsv")),
      meta_path      = file.path(RESULTS, paste0("meta_analysis__", db, ".tsv")),
      row_filter_fn  = cfg$filter_fn,
      concordance_q  = thresholds$star1
    )
    n_rows <- nrow(dat$coef)
    if (n_rows == 0) { message("  skipped: 0 rows after filter"); next }
    height <- min(40, max(8, 0.18 * n_rows + 3))

    for (order_mode in c("padj", "clust")) {
      ht <- build_heatmap(
        dat,
        order_mode     = order_mode,
        title_prefix   = sprintf("%s (%d genera)", toupper(db), n_rows),
        row_label_fn   = strip_genus_prefix,
        thresholds     = thresholds,
        meta_indicator = cfg$indicator
      )
      base <- file.path(HEATDIR,
                        sprintf("heatmap_%s%s_by_%s", db, cfg$tag, order_mode))
      draw_to_device <- function() {
        ComplexHeatmap::draw(
          ht,
          merge_legend           = TRUE,
          heatmap_legend_side    = "right",
          annotation_legend_side = "right"
        )
      }

      # Vector PDF (existing main output)
      out_pdf <- paste0(base, ".pdf")
      pdf(out_pdf, width = 11, height = height)
      draw_to_device()
      dev.off()
      message("  wrote ", out_pdf, "  (", n_rows, " rows)")

      # 300 DPI PNG for review / draft inclusion
      out_png <- paste0(base, ".png")
      png(out_png, width = 11, height = height, units = "in", res = 300,
          type = "cairo")
      draw_to_device()
      dev.off()
      message("  wrote ", out_png)

      # 600 DPI LZW-compressed TIFF for journal submission
      out_tiff <- paste0(base, ".tiff")
      tiff(out_tiff, width = 11, height = height, units = "in", res = 600,
           compression = "lzw", type = "cairo")
      draw_to_device()
      dev.off()
      message("  wrote ", out_tiff)
    }
  }
}

message("\nDone.")
message("Main figure (directional consistency): heatmap_{silva,ehomd}_v2_by_{padj,clust}.{pdf,png,tiff}")
message("Supplementary (all meta features):     heatmap_{silva,ehomd}_by_{padj,clust}.{pdf,png,tiff}")
