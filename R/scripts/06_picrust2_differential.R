# 06_picrust2_differential.R — Per-study MaAsLin2 on PICRUSt2 functional predictions +
# random-effects meta-analysis per disease. Mirrors scripts 04 + 05 but on functional
# features (MetaCyc pathways and KEGG KOs) instead of taxa.
#
# Inputs (per study):
#   {DATA_ROOT}/{dir}/exported-picrust2/pathway_abundance.tsv   # MetaCyc pathway × samples
#   {DATA_ROOT}/{dir}/exported-picrust2/ko_metagenome.tsv       # KEGG KO × samples
#   cache/ps_list.rds                                            # for group metadata
# Outputs:
#   results/maaslin2_picrust2/<study>__{pathway,ko}/all_results.tsv
#   results/picrust2_all_studies__{pathway,ko}.tsv   (pooled per-study MaAsLin2)
#   results/meta_analysis_picrust2__{pathway,ko}.tsv  (metafor REML per disease)
#
# Notes
# - Feature types are analyzed separately (not per DB — PICRUSt2 is DB-agnostic).
# - Metadata (group assignment) is taken from the SILVA phyloseq object so samples
#   match the taxonomic analyses in scripts 02–05.
# - TSS + LOG normalization is used (same as script 04) — PICRUSt2 abundances are
#   predicted counts, so relative-abundance + log transform is reasonable.

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))

ps_list <- readRDS(file.path(CACHE, "ps_list.rds"))

# Read a PICRUSt2 biom-exported TSV (header: "# Constructed from biom file" then "#OTU ID").
read_picrust2 <- function(path) {
  df <- readr::read_tsv(path, skip = 1, show_col_types = FALSE)
  colnames(df)[1] <- "feature"
  df
}

feature_files <- list(
  pathway = "pathway_abundance.tsv",
  ko      = "ko_metagenome.tsv"
)

# MaAsLin2 per study x feature type
for (s in studies$study) {
  study_dir <- file.path(DATA_ROOT, studies$dir[studies$study == s])

  ps <- ps_list[[paste0(s, "__silva")]]
  if (is.null(ps)) {
    message("No phyloseq for ", s, "; skip.")
    next
  }
  meta_raw <- as(phyloseq::sample_data(ps), "data.frame")
  grp_col  <- grep("^group$", colnames(meta_raw), ignore.case = TRUE, value = TRUE)[1]
  meta_df  <- data.frame(
    group     = as.character(meta_raw[[grp_col]]),
    row.names = phyloseq::sample_names(ps),
    stringsAsFactors = FALSE
  )

  for (ftype in names(feature_files)) {
    path <- file.path(study_dir, "exported-picrust2", feature_files[[ftype]])
    if (!file.exists(path)) {
      message("Skip (missing ", ftype, "): ", s)
      next
    }
    message("MaAsLin2 PICRUSt2: ", s, " [", ftype, "]")

    mat <- read_picrust2(path)
    feat_mat <- as.data.frame(mat[, -1], stringsAsFactors = FALSE)
    rownames(feat_mat) <- mat$feature
    # MaAsLin2 expects rows = samples, cols = features
    input_data <- as.data.frame(t(feat_mat), stringsAsFactors = FALSE)

    keep <- intersect(rownames(input_data), rownames(meta_df))
    if (length(keep) < 4) {
      message("  <4 matching samples (", length(keep), "); skip.")
      next
    }
    input_data_s <- input_data[keep, , drop = FALSE]
    meta_df_s    <- meta_df[keep, , drop = FALSE]

    out_dir <- file.path(RESULTS, "maaslin2_picrust2", paste0(s, "__", ftype))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    Maaslin2::Maaslin2(
      input_data     = input_data_s,
      input_metadata = meta_df_s,
      output         = out_dir,
      fixed_effects  = "group",
      normalization  = "TSS",
      transform      = "LOG",
      min_prevalence = 0.1,
      min_abundance  = 0.0001,   # drop features present but vanishingly rare
      plot_heatmap   = TRUE,
      plot_scatter   = FALSE
    )

    # Back-map MaAsLin2's sanitized feature names (e.g. X1CMET2.PWY) to the
    # original PICRUSt2 IDs (1CMET2-PWY) so results stay human-readable and
    # searchable against MetaCyc / KEGG.
    id_map <- tibble::tibble(
      original  = mat$feature,
      sanitized = make.names(mat$feature)
    )
    for (f in c("all_results.tsv", "significant_results.tsv")) {
      fp <- file.path(out_dir, f)
      if (!file.exists(fp)) next
      res <- readr::read_tsv(fp, show_col_types = FALSE)
      res <- dplyr::left_join(res, id_map, by = c("feature" = "sanitized")) |>
        dplyr::mutate(feature = dplyr::coalesce(original, feature)) |>
        dplyr::select(-original)
      readr::write_tsv(res, fp)
    }
  }
}

# Meta-analysis per feature type x disease
for (ftype in names(feature_files)) {
  message("\n=== PICRUSt2 meta-analysis: ", ftype, " ===")

  maaslin_files <- tibble::tibble(
    study   = studies$study,
    disease = studies$disease,
    path    = file.path(RESULTS, "maaslin2_picrust2",
                        paste0(studies$study, "__", ftype),
                        "all_results.tsv")
  ) |> dplyr::filter(file.exists(path))

  if (nrow(maaslin_files) == 0) {
    message("  No MaAsLin2 outputs found for ", ftype, "; skip.")
    next
  }

  all_da <- purrr::pmap_dfr(maaslin_files, function(study, disease, path) {
    readr::read_tsv(path, show_col_types = FALSE) |>
      dplyr::mutate(study = study, disease = disease)
  })

  readr::write_tsv(all_da,
                   file.path(RESULTS, paste0("picrust2_all_studies__", ftype, ".tsv")))

  meta_results <- all_da |>
    dplyr::group_by(disease, feature) |>
    dplyr::filter(dplyr::n() >= 2) |>
    dplyr::group_modify(~{
      fit <- tryCatch(
        metafor::rma(yi = .x$coef, sei = .x$stderr, method = "DL"),
        error = function(e) NULL
      )
      if (is.null(fit)) return(tibble::tibble())
      tibble::tibble(k = fit$k, estimate = fit$beta[1], se = fit$se,
                     p = fit$pval, I2 = fit$I2)
    }) |>
    dplyr::ungroup() |>
    dplyr::group_by(disease) |>
    dplyr::mutate(padj = p.adjust(p, "BH")) |>   # BH within disease
    dplyr::ungroup() |>
    dplyr::arrange(padj)

  readr::write_tsv(meta_results,
                   file.path(RESULTS, paste0("meta_analysis_picrust2__", ftype, ".tsv")))
  message("  ", nrow(meta_results), " feature-disease pairs tested. ",
          "Top results written to meta_analysis_picrust2__", ftype, ".tsv")
}

message("\nPICRUSt2 differential abundance done. ",
        "See results/maaslin2_picrust2/ and results/meta_analysis_picrust2__{pathway,ko}.tsv")
