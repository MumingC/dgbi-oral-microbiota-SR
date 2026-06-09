# 05b_sensitivity_amplicon.R — Sensitivity analysis: meta-analysis restricted to
# V3-V4 studies only, to assess whether the main meta-analysis findings are driven
# by amplicon-region heterogeneity (Kawar 2021 = V1-V3; Tang 2023 = V4-V5).
#
# Strategy:
#   - Drop Kawar_2021 (V1-V3) and Tang_2023 (V4-V5) from the pooled per-study table.
#   - Re-run random-effects meta-analysis on V3-V4 studies only.
#   - Compare results to the main analysis (05_meta_analysis.R).
#
# Consequence for IBS-D:
#   - IBS-D V3-V4 studies: Li_2025 only (k=1). Meta-analysis requires k>=2.
#   - Therefore IBS-D taxa sensitivity result = "cannot assess; only one V3-V4 study."
#   - This is itself a finding: the IBS-D taxa meta signal rests entirely on the
#     pairing of Li (V3-V4) and Tang (V4-V5) — cross-region.
#
# Consequence for GERD:
#   - GERD V3-V4 studies: Ziganshina_2020, Qian_2023, Hao_2022, Zheng_2024 (k=4).
#   - Kawar_2021 removed. Meta-analysis feasible and comparable to main.
#
# Outputs:
#   results/sensitivity_amplicon/meta_analysis__{silva,ehomd}_v3v4only.tsv
#   results/sensitivity_amplicon/meta_analysis_picrust2__{pathway,ko}_v3v4only.tsv

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))

SENS_DIR <- file.path(RESULTS, "sensitivity_amplicon")
dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)

V3V4_STUDIES <- studies$study[studies$region == "V3-V4"]
EXCL <- setdiff(studies$study, V3V4_STUDIES)
message("Sensitivity analysis: V3-V4 only")
message("  Included: ", paste(V3V4_STUDIES, collapse = ", "))
message("  Excluded: ", paste(EXCL, collapse = ", "))

run_sensitivity_meta <- function(pooled_path, out_path, label) {
  if (!file.exists(pooled_path)) {
    message("  Missing: ", pooled_path, " — skip")
    return(invisible(NULL))
  }
  all_da <- readr::read_tsv(pooled_path, show_col_types = FALSE) |>
    dplyr::filter(study %in% V3V4_STUDIES)

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

  readr::write_tsv(meta_results, out_path)

  n_sig <- sum(meta_results$padj < 0.05, na.rm = TRUE)
  by_dis <- meta_results |>
    dplyr::group_by(disease) |>
    dplyr::summarise(n_tested = dplyr::n(),
                     n_padj05 = sum(padj < 0.05, na.rm = TRUE),
                     .groups = "drop")
  message("\n  === ", label, " ===")
  message("  Total tested: ", nrow(meta_results), " | padj<0.05: ", n_sig)
  print(as.data.frame(by_dis))
  invisible(meta_results)
}

# Taxa
for (db in c("silva", "ehomd")) {
  run_sensitivity_meta(
    pooled_path = file.path(RESULTS, paste0("maaslin2_all_studies__", db, ".tsv")),
    out_path    = file.path(SENS_DIR, paste0("meta_analysis__", db, "_v3v4only.tsv")),
    label       = paste("Taxa", toupper(db), "V3-V4 only")
  )
}

# PICRUSt2
for (ftype in c("pathway", "ko")) {
  run_sensitivity_meta(
    pooled_path = file.path(RESULTS, paste0("picrust2_all_studies__", ftype, ".tsv")),
    out_path    = file.path(SENS_DIR, paste0("meta_analysis_picrust2__", ftype, "_v3v4only.tsv")),
    label       = paste("PICRUSt2", ftype, "V3-V4 only")
  )
}

message("\nSensitivity analysis done. See results/sensitivity_amplicon/")
message("Note: IBS-D will have k<2 for taxa — this is expected and itself a finding.")
message("Compare these results against the main meta_analysis_*.tsv outputs.")
