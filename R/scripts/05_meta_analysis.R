# 05_meta_analysis.R — Random-effects meta-analysis of per-genus effects within each disease.
# Inputs: per-study MaAsLin2 all_results.tsv for BOTH SILVA and eHOMD.
# Pools log2FC estimates via metafor::rma() separately per database.
#
# BH correction is applied WITHIN disease (not globally) so that GERD and IBS-D
# multiple-testing burdens are kept separate. This is more interpretable and
# standard for disease-stratified meta-analyses.

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))

for (db in c("silva", "ehomd")) {
  message("\n=== Meta-analysis: ", db, " ===")

  maaslin_files <- tibble::tibble(
    study   = studies$study,
    disease = studies$disease,
    path    = file.path(RESULTS, "maaslin2",
                        paste0(studies$study, "__", db),
                        "all_results.tsv")
  ) |> dplyr::filter(file.exists(path))

  if (nrow(maaslin_files) == 0) {
    message("  No MaAsLin2 outputs found for ", db, "; skip.")
    next
  }

  all_da <- purrr::pmap_dfr(maaslin_files, function(study, disease, path) {
    readr::read_tsv(path, show_col_types = FALSE) |>
      dplyr::mutate(study = study, disease = disease)
  })

  readr::write_tsv(all_da,
                   file.path(RESULTS, paste0("maaslin2_all_studies__", db, ".tsv")))

  # Meta-analysis per disease × genus (keep genera present in >=2 studies)
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
                   file.path(RESULTS, paste0("meta_analysis__", db, ".tsv")))
  message("  ", nrow(meta_results), " genus-disease pairs tested. ",
          "Top results written to meta_analysis__", db, ".tsv")
}

message("\nMeta-analysis done. See results/meta_analysis__{silva,ehomd}.tsv")
