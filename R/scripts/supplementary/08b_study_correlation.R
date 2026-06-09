# 08b_study_correlation.R — Pairwise per-study coefficient correlation.
#
# Sanity check for "is Zheng 2024 an outlier in the Upper-GI group?" Builds a
# 7x7 Spearman correlation matrix of per-study MaAsLin2 coef vectors, paired
# by shared genera (pairwise complete observations). A study that correlates
# strongly with its disease-siblings and weakly with out-of-disease studies
# supports the current grouping; a study with low correlations across the
# board supports teasing it out or dropping it.
#
# Outputs (figures/heatmap/):
#   study_coef_correlation_silva.pdf    # vector, main check
#   study_coef_correlation_silva.png    # 300 DPI, for review / draft
#   study_coef_correlation_silva.tiff   # 600 DPI LZW, journal submission
#   study_coef_correlation_ehomd.{pdf,png,tiff}   # cross-database confirmation
#
# Each cell:  Spearman r  (upper-right half)
#             N overlap   (lower-left half)
# Diagonal shows self-comparison (r=1, N=#features in that study).

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
}
if (!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

read_safe <- function(path) {
  raw <- readBin(path, "raw", file.info(path)$size)
  raw <- raw[raw != as.raw(0)]
  readr::read_tsv(I(rawToChar(raw)), show_col_types = FALSE)
}

STUDY_ORDER <- c("Ziganshina_2020", "Qian_2023", "Hao_2022",
                 "Kawar_2021", "Zheng_2024",
                 "Li_2025", "Tang_2023")
STUDY_DISEASE <- setNames(
  c("GERD", "GERD", "GERD", "GERD", "LPRD", "IBS-D", "IBS-D"),
  STUDY_ORDER
)

HEATDIR <- file.path(FIGURES, "heatmap")
dir.create(HEATDIR, recursive = TRUE, showWarnings = FALSE)

build_corr <- function(db) {
  per_study <- read_safe(
    file.path(RESULTS, paste0("maaslin2_all_studies__", db, ".tsv"))
  )
  wide <- per_study %>%
    dplyr::select(feature, study, coef) %>%
    tidyr::pivot_wider(names_from = study, values_from = coef) %>%
    tibble::column_to_rownames("feature") %>%
    as.matrix()
  wide <- wide[, intersect(STUDY_ORDER, colnames(wide)), drop = FALSE]

  n_studies <- ncol(wide)
  r_mat <- matrix(NA_real_, n_studies, n_studies,
                  dimnames = list(colnames(wide), colnames(wide)))
  n_mat <- matrix(NA_integer_, n_studies, n_studies,
                  dimnames = list(colnames(wide), colnames(wide)))

  for (i in seq_len(n_studies)) {
    for (j in seq_len(n_studies)) {
      xi <- wide[, i]; xj <- wide[, j]
      ok <- !is.na(xi) & !is.na(xj)
      n_mat[i, j] <- sum(ok)
      if (sum(ok) >= 3) {
        r_mat[i, j] <- suppressWarnings(
          cor(xi[ok], xj[ok], method = "spearman")
        )
      }
    }
  }
  list(r = r_mat, n = n_mat)
}

draw_corr <- function(cc, db, out_pdf, out_png = NULL, out_tiff = NULL) {
  r_mat <- cc$r
  n_mat <- cc$n
  col_fun <- circlize::colorRamp2(
    c(-1, 0, 1),
    c("#3970b0", "white", "#d4652a")
  )
  disease_vec <- STUDY_DISEASE[rownames(r_mat)]
  row_ann <- ComplexHeatmap::rowAnnotation(
    disease = disease_vec,
    col = list(disease = c(GERD = "#f4a261", LPRD = "#e76f51",
                           `IBS-D` = "#2a9d8f")),
    show_annotation_name = FALSE
  )
  col_ann <- ComplexHeatmap::HeatmapAnnotation(
    disease = disease_vec,
    col = list(disease = c(GERD = "#f4a261", LPRD = "#e76f51",
                           `IBS-D` = "#2a9d8f")),
    show_annotation_name = FALSE
  )
  # layer_fun: upper-right = r, lower-left = N, diagonal = both
  layer_fun <- function(j, i, x, y, width, height, fill) {
    r_vec <- ComplexHeatmap::pindex(r_mat, i, j)
    n_vec <- ComplexHeatmap::pindex(n_mat, i, j)
    upper <- j > i
    lower <- j < i
    diag  <- j == i
    # upper: r
    if (any(upper)) {
      lab <- sprintf("%.2f", r_vec[upper])
      lab[is.na(r_vec[upper])] <- ""
      grid::grid.text(lab, x = x[upper], y = y[upper],
                      gp = grid::gpar(fontsize = 8))
    }
    # lower: N
    if (any(lower)) {
      grid::grid.text(n_vec[lower], x = x[lower], y = y[lower],
                      gp = grid::gpar(fontsize = 8, col = "grey30"))
    }
    # diagonal: "n=N"
    if (any(diag)) {
      grid::grid.text(sprintf("n=%d", n_vec[diag]),
                      x = x[diag], y = y[diag],
                      gp = grid::gpar(fontsize = 8, fontface = "bold"))
    }
  }
  ht <- ComplexHeatmap::Heatmap(
    r_mat,
    name             = "Spearman r",
    col              = col_fun,
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    top_annotation   = col_ann,
    left_annotation  = row_ann,
    row_names_gp     = grid::gpar(fontsize = 9),
    column_names_gp  = grid::gpar(fontsize = 9),
    column_names_rot = 45,
    layer_fun        = layer_fun,
    rect_gp          = grid::gpar(col = "white", lwd = 0.5),
    heatmap_legend_param = list(
      at = c(-1, -0.5, 0, 0.5, 1), direction = "vertical"
    ),
    column_title = sprintf(
      "%s \u2014 per-study coef correlation (upper: r, lower: N overlap)",
      toupper(db)
    ),
    column_title_gp = grid::gpar(fontsize = 11, fontface = "bold")
  )
  draw_to_device <- function() {
    ComplexHeatmap::draw(ht, merge_legend = TRUE,
                         heatmap_legend_side = "right",
                         annotation_legend_side = "right")
  }

  # Vector PDF (existing main output)
  pdf(out_pdf, width = 9, height = 7)
  draw_to_device()
  dev.off()
  message("wrote ", out_pdf)

  # 300 DPI PNG for review / draft inclusion
  if (!is.null(out_png)) {
    png(out_png, width = 9, height = 7, units = "in", res = 300,
        type = "cairo")
    draw_to_device()
    dev.off()
    message("wrote ", out_png)
  }

  # 600 DPI LZW-compressed TIFF for journal submission
  if (!is.null(out_tiff)) {
    tiff(out_tiff, width = 9, height = 7, units = "in", res = 600,
         compression = "lzw", type = "cairo")
    draw_to_device()
    dev.off()
    message("wrote ", out_tiff)
  }
}

for (db in c("silva", "ehomd")) {
  message("\n=== ", db, " ===")
  cc <- build_corr(db)
  # Print matrix to console for quick inspection.
  message("Spearman r:")
  print(round(cc$r, 2))
  message("N overlap:")
  print(cc$n)

  out_pdf  <- file.path(HEATDIR, sprintf("study_coef_correlation_%s.pdf",  db))
  out_png  <- file.path(HEATDIR, sprintf("study_coef_correlation_%s.png",  db))
  out_tiff <- file.path(HEATDIR, sprintf("study_coef_correlation_%s.tiff", db))
  draw_corr(cc, db, out_pdf, out_png, out_tiff)

  # Also write tidy TSVs for the record.
  r_df <- as.data.frame(cc$r) %>% tibble::rownames_to_column("study")
  n_df <- as.data.frame(cc$n) %>% tibble::rownames_to_column("study")
  readr::write_tsv(r_df,
    file.path(RESULTS, sprintf("study_coef_correlation_%s.tsv", db)))
  readr::write_tsv(n_df,
    file.path(RESULTS, sprintf("study_coef_overlap_%s.tsv", db)))
}

message("\nDone. See figures/heatmap/study_coef_correlation_{silva,ehomd}.{pdf,png,tiff}")
