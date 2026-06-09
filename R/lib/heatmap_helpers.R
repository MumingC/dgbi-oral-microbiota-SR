# heatmap_helpers.R — Shared code for Duvallet-style effect-size heatmaps
# (taxa in scripts/08_heatmap.R, PICRUSt2 pathways in scripts/09_heatmap_picrust2.R).
#
# Exposed API:
#   install_heatmap_deps()                       — lazy-install ComplexHeatmap + circlize
#   read_safe(path)                              — NUL-safe TSV reader
#   panels                                       — named list of panel definitions
#   STUDY_DISEASE, DISEASE_COLORS                — per-study colour annotation
#   build_heatmap_data(per_study_path, meta_path, row_filter_fn)
#   build_heatmap(dat, order_mode, title_prefix, row_label_fn, thresholds)
#
# "dat" is a list with:
#   coef, qval, col_panel, col_is_meta, col_disease, best_padj

install_heatmap_deps <- function() {
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
}

# Strip embedded NULs before parsing (some TSVs have trailing NUL padding).
read_safe <- function(path) {
  raw <- readBin(path, "raw", file.info(path)$size)
  raw <- raw[raw != as.raw(0)]
  readr::read_tsv(I(rawToChar(raw)), show_col_types = FALSE)
}

# ---- Panel + colour constants ----------------------------------------------
panels <- list(
  "Upper-GI" = list(
    studies      = c("Ziganshina_2020", "Qian_2023", "Hao_2022",
                     "Kawar_2021", "Zheng_2024"),
    meta_disease = "GERD",
    meta_label   = "GERD Meta"
  ),
  "IBS-D" = list(
    studies      = c("Li_2025", "Tang_2023"),
    meta_disease = "IBS-D",
    meta_label   = "IBS-D Meta"
  )
)

STUDY_DISEASE <- c(
  Ziganshina_2020 = "GERD", Qian_2023 = "GERD",
  Hao_2022        = "GERD", Kawar_2021 = "GERD",
  Zheng_2024      = "LPRD",
  Li_2025         = "IBS-D", Tang_2023 = "IBS-D"
)
DISEASE_COLORS <- c(GERD = "#f4a261", LPRD = "#e76f51", `IBS-D` = "#2a9d8f")

# ---- Directional concordance ------------------------------------------------
# For one genus, given per-study coefs and qvals restricted to a panel's
# meta_disease studies, decide if there is directional concordance:
#   >=2 studies share the same sign of coef, AND
#   >=1 of those agreeing studies has qval < q_thresh.
# Returns list with: pass (logical), direction ("up"/"dn"/"mixed"/"na").
concordance_one <- function(coef_vec, qval_vec, q_thresh = 0.25) {
  valid <- !is.na(coef_vec)
  cv <- coef_vec[valid]; qv <- qval_vec[valid]
  if (length(cv) < 2) return(list(pass = FALSE, direction = "na"))
  up_idx <- which(cv > 0); dn_idx <- which(cv < 0)
  up_pass <- length(up_idx) >= 2 && any(!is.na(qv[up_idx]) & qv[up_idx] < q_thresh)
  dn_pass <- length(dn_idx) >= 2 && any(!is.na(qv[dn_idx]) & qv[dn_idx] < q_thresh)
  if (up_pass && !dn_pass) return(list(pass = TRUE,  direction = "up"))
  if (dn_pass && !up_pass) return(list(pass = TRUE,  direction = "dn"))
  if (up_pass && dn_pass)  return(list(pass = FALSE, direction = "mixed"))
  list(pass = FALSE, direction = "na")
}

# Per-panel concordance call for every feature. Restricted to studies whose
# disease matches the panel's meta_disease, so LPRD (Zheng_2024) does NOT
# contribute to GERD concordance.
# Returns list: panel -> named character vector (feature -> direction).
compute_panel_concordance <- function(per_study_df, q_thresh = 0.25) {
  out <- list()
  for (pn in names(panels)) {
    dis      <- panels[[pn]]$meta_disease
    dstudies <- names(STUDY_DISEASE)[STUDY_DISEASE == dis]
    dstudies <- intersect(dstudies, unique(per_study_df$study))
    sub <- per_study_df[per_study_df$study %in% dstudies,
                        c("feature", "study", "coef", "qval"), drop = FALSE]
    if (nrow(sub) == 0) { out[[pn]] <- character(); next }
    sub_split <- split(sub, sub$feature)
    dv <- vapply(sub_split, function(df) {
      concordance_one(df$coef, df$qval, q_thresh = q_thresh)$direction
    }, character(1))
    out[[pn]] <- dv
  }
  out
}

# Row filter: keep genera that are concordant in >=1 panel.
row_filter_concordance <- function(meta_df, per_study_df, q_thresh = 0.25) {
  panel_dir <- compute_panel_concordance(per_study_df, q_thresh = q_thresh)
  feats <- character()
  for (pn in names(panel_dir)) {
    dv <- panel_dir[[pn]]
    feats <- c(feats, names(dv)[dv %in% c("up", "dn")])
  }
  unique(feats)
}

# ---- Build the full data object from per-study + meta TSVs -----------------
# row_filter_fn(meta_df, per_study_df) returns a character vector of feature
# names to keep. Filters that don't need per_study_df may ignore the second arg.
build_heatmap_data <- function(per_study_path, meta_path, row_filter_fn,
                               concordance_q = 0.25) {
  per_study <- read_safe(per_study_path)
  meta      <- read_safe(meta_path)

  row_feats <- row_filter_fn(meta, per_study)
  if (length(row_feats) == 0) stop("row_filter_fn returned 0 features.")

  long <- per_study %>%
    dplyr::filter(feature %in% row_feats) %>%
    dplyr::select(feature, study, coef, qval)

  mat_coef <- long %>%
    dplyr::select(feature, study, coef) %>%
    tidyr::pivot_wider(names_from = study, values_from = coef) %>%
    tibble::column_to_rownames("feature") %>%
    as.matrix()
  mat_qval <- long %>%
    dplyr::select(feature, study, qval) %>%
    tidyr::pivot_wider(names_from = study, values_from = qval) %>%
    tibble::column_to_rownames("feature") %>%
    as.matrix()

  study_order <- unlist(lapply(panels, `[[`, "studies"), use.names = FALSE)
  study_order <- intersect(study_order, colnames(mat_coef))
  mat_coef <- mat_coef[, study_order, drop = FALSE]
  mat_qval <- mat_qval[, study_order, drop = FALSE]

  build_meta_col <- function(dis) {
    m <- meta %>% dplyr::filter(disease == dis)
    vec  <- setNames(m$estimate, m$feature)
    padj <- setNames(m$padj,     m$feature)
    list(val = vec[rownames(mat_coef)], padj = padj[rownames(mat_coef)])
  }
  meta_cols <- list()
  for (pn in names(panels)) {
    d   <- panels[[pn]]$meta_disease
    lbl <- panels[[pn]]$meta_label
    if (any(meta$disease == d)) meta_cols[[lbl]] <- build_meta_col(d)
  }

  full_coef <- NULL; full_qval <- NULL
  col_panel <- c(); col_is_meta <- c(); col_disease <- c()
  for (pn in names(panels)) {
    pstuds <- intersect(panels[[pn]]$studies, colnames(mat_coef))
    if (length(pstuds)) {
      full_coef <- cbind(full_coef, mat_coef[, pstuds, drop = FALSE])
      full_qval <- cbind(full_qval, mat_qval[, pstuds, drop = FALSE])
      col_panel   <- c(col_panel,   rep(pn, length(pstuds)))
      col_is_meta <- c(col_is_meta, rep(FALSE, length(pstuds)))
      col_disease <- c(col_disease, STUDY_DISEASE[pstuds])
    }
    mlbl <- panels[[pn]]$meta_label
    if (!is.null(meta_cols[[mlbl]])) {
      full_coef <- cbind(full_coef, meta_cols[[mlbl]]$val)
      full_qval <- cbind(full_qval, meta_cols[[mlbl]]$padj)
      colnames(full_coef)[ncol(full_coef)] <- mlbl
      colnames(full_qval)[ncol(full_qval)] <- mlbl
      col_panel   <- c(col_panel,   pn)
      col_is_meta <- c(col_is_meta, TRUE)
      col_disease <- c(col_disease, panels[[pn]]$meta_disease)
    }
  }

  meta_wide <- meta %>%
    dplyr::select(feature, disease, padj) %>%
    tidyr::pivot_wider(names_from = disease, values_from = padj,
                       names_prefix = "padj_") %>%
    tibble::column_to_rownames("feature")
  best_padj <- apply(meta_wide[rownames(full_coef), , drop = FALSE], 1,
                     function(x) suppressWarnings(min(x, na.rm = TRUE)))
  best_padj[is.infinite(best_padj)] <- NA_real_

  panel_concordance <- compute_panel_concordance(per_study,
                                                 q_thresh = concordance_q)

  list(
    coef              = full_coef,
    qval              = full_qval,
    col_panel         = col_panel,
    col_is_meta       = col_is_meta,
    col_disease       = col_disease,
    best_padj         = best_padj,
    panel_concordance = panel_concordance
  )
}

# ---- Build a ComplexHeatmap object ------------------------------------------
# row_label_fn(rownames) returns the display labels. Defaults to identity.
# thresholds is a list with: star1, star2, meta_star1, meta_star2, meta_sig, color_limit
#   star1 / star2          — per-study cell thresholds (* / **)
#   meta_star1 / meta_star2 — Meta cell thresholds (* / **). Set meta_star1 = NA
#                            to suppress the * tier and show only ** at meta padj
#                            < meta_star2 — the stricter scheme appropriate for a
#                            meta-analytic inference column.
build_heatmap <- function(dat,
                          order_mode,
                          title_prefix,
                          row_label_fn = identity,
                          thresholds = list(star1 = 0.25, star2 = 0.05,
                                            meta_star1 = NA, meta_star2 = 0.05,
                                            meta_sig = 0.05, color_limit = 3),
                          meta_indicator = c("border", "triangle", "none", "star")) {
  meta_indicator <- match.arg(meta_indicator)
  coef_mat <- dat$coef
  qval_mat <- dat$qval

  if (order_mode == "padj") {
    ord <- order(dat$best_padj, na.last = TRUE)
    coef_mat <- coef_mat[ord, , drop = FALSE]
    qval_mat <- qval_mat[ord, , drop = FALSE]
    cluster_rows <- FALSE
    title_suffix <- "ordered by min meta padj"
  } else if (order_mode == "clust") {
    mat_for_clust <- coef_mat
    mat_for_clust[is.na(mat_for_clust)] <- 0
    cluster_rows <- hclust(dist(mat_for_clust), method = "average")
    title_suffix <- "hierarchical clustering of coef profile"
  } else {
    stop("unknown order_mode: ", order_mode)
  }

  CLIM <- thresholds$color_limit
  col_fun <- circlize::colorRamp2(
    c(-CLIM, 0, CLIM),
    c("#3970b0", "white", "#d4652a")
  )

  s1   <- thresholds$star1
  s2   <- thresholds$star2
  ms1  <- thresholds$meta_star1   # NA => suppress * tier on Meta cells
  ms2  <- thresholds$meta_star2
  msig <- thresholds$meta_sig
  if (is.null(ms1)) ms1 <- NA_real_
  if (is.null(ms2)) ms2 <- 0.05
  panel_concordance <- dat$panel_concordance
  layer_fun <- function(j, i, x, y, width, height, fill) {
    is_meta_vec <- dat$col_is_meta[j]
    q_vec       <- ComplexHeatmap::pindex(qval_mat, i, j)
    c_vec       <- ComplexHeatmap::pindex(coef_mat, i, j)

    # Per-study asterisks (unchanged).
    s2_mask <- !is_meta_vec & !is.na(q_vec) & q_vec < s2
    if (any(s2_mask)) grid::grid.text("**", x = x[s2_mask], y = y[s2_mask],
                                      gp = grid::gpar(fontsize = 7))
    s1_mask <- !is_meta_vec & !is.na(q_vec) & q_vec >= s2 & q_vec < s1
    if (any(s1_mask)) grid::grid.text("*", x = x[s1_mask], y = y[s1_mask],
                                      gp = grid::gpar(fontsize = 7))

    # Meta-cell indicator.
    if (meta_indicator == "border") {
      mbox <- is_meta_vec & !is.na(q_vec) & q_vec < msig
      if (any(mbox)) grid::grid.rect(
        x = x[mbox], y = y[mbox],
        width = width[mbox], height = height[mbox],
        gp = grid::gpar(col = "black", lwd = 1.8, fill = NA)
      )
    } else if (meta_indicator == "triangle") {
      mm <- which(is_meta_vec)
      if (length(mm)) {
        pv <- dat$col_panel[j[mm]]
        fv <- rownames(coef_mat)[i[mm]]
        mc <- c_vec[mm]
        dir_char <- vapply(seq_along(pv), function(k) {
          pc <- panel_concordance[[pv[k]]]
          if (!is.null(pc) && fv[k] %in% names(pc)) pc[[fv[k]]] else "na"
        }, character(1))
        up_ok <- dir_char == "up" & !is.na(mc) & mc > 0
        dn_ok <- dir_char == "dn" & !is.na(mc) & mc < 0
        if (any(up_ok)) grid::grid.points(
          x[mm][up_ok], y[mm][up_ok],
          pch = 24, size = grid::unit(2.5, "mm"),
          gp = grid::gpar(fill = "black", col = "black")
        )
        if (any(dn_ok)) grid::grid.points(
          x[mm][dn_ok], y[mm][dn_ok],
          pch = 25, size = grid::unit(2.5, "mm"),
          gp = grid::gpar(fill = "black", col = "black")
        )
      }
    } else if (meta_indicator == "star") {
      # Pooled-column significance. Uses meta-specific thresholds (ms1, ms2)
      # which are typically stricter than per-study (s1, s2). Default config:
      # ms1 = NA (no * tier), ms2 = 0.05 — i.e. ** only at meta padj < 0.05.
      ms2_mask <- is_meta_vec & !is.na(q_vec) & q_vec < ms2
      if (any(ms2_mask)) grid::grid.text("**", x = x[ms2_mask], y = y[ms2_mask],
                                         gp = grid::gpar(fontsize = 7))
      if (!is.na(ms1)) {
        ms1_mask <- is_meta_vec & !is.na(q_vec) & q_vec >= ms2 & q_vec < ms1
        if (any(ms1_mask)) grid::grid.text("*", x = x[ms1_mask], y = y[ms1_mask],
                                           gp = grid::gpar(fontsize = 7))
      }
    }
    # meta_indicator == "none": nothing.
  }

  top_ann <- ComplexHeatmap::HeatmapAnnotation(
    disease = dat$col_disease,
    col     = list(disease = DISEASE_COLORS),
    show_annotation_name = FALSE,
    annotation_legend_param = list(disease = list(title = "Disease"))
  )

  ComplexHeatmap::Heatmap(
    coef_mat,
    name             = "MaAsLin2 coef",
    col              = col_fun,
    cluster_rows     = cluster_rows,
    cluster_columns  = FALSE,
    column_split     = factor(dat$col_panel, levels = names(panels)),
    column_gap       = grid::unit(3, "mm"),
    column_title_gp  = grid::gpar(fontsize = 12, fontface = "bold"),
    top_annotation   = top_ann,
    row_labels       = row_label_fn(rownames(coef_mat)),
    row_names_gp     = grid::gpar(fontsize = 7),
    row_names_side   = "left",
    # Widen row-label allocation so long MetaCyc English names (e.g. "superpathway
    # of geranylgeranyl diphosphate biosynthesis I (via mevalonate)") don't clip.
    row_names_max_width = grid::unit(14, "cm"),
    column_names_gp  = grid::gpar(fontsize = 8),
    column_names_rot = 45,
    layer_fun        = layer_fun,
    na_col           = "grey92",
    rect_gp          = grid::gpar(col = "white", lwd = 0.3),
    border           = TRUE,
    heatmap_legend_param = list(
      title     = "MaAsLin2 coef\n(+ = up in disease)",
      at        = c(-CLIM, 0, CLIM),
      direction = "vertical"
    ),
    column_title = sprintf("%s \u2014 %s", title_prefix, title_suffix)
  )
}
