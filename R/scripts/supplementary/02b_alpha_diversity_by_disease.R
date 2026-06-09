# 02b_alpha_diversity_by_disease.R
# Regenerates figures/alpha_diversity.pdf with studies ordered and grouped by disease.
# Replaces the Apr-15 version which was not organised by disease.
# Run from dgbi_exports/R/: Rscript scripts/02b_alpha_diversity_by_disease.R

source(here::here("00_setup.R"))

RESULTS <- here::here("results")
FIGURES <- here::here("figures")

alpha_df <- readr::read_tsv(file.path(RESULTS, "alpha_diversity.tsv"),
                            show_col_types = FALSE)

alpha_df$case <- ifelse(
  tolower(alpha_df$group) %in% c("control", "healthy", "hc"),
  "Control", "Case"
)

# ---- Order studies by disease: GERD → LPRD → IBS-D --------------------------
disease_order <- c("GERD", "LPRD", "IBS-D")
study_order <- alpha_df |>
  dplyr::distinct(study, disease) |>
  dplyr::mutate(disease = factor(disease, levels = disease_order)) |>
  dplyr::arrange(disease) |>
  dplyr::pull(study)

alpha_df <- alpha_df |>
  dplyr::mutate(
    study   = factor(study, levels = study_order),
    disease = factor(disease, levels = disease_order),
    # Study label for facet strip: "Author Year\n(DISEASE)"
    study_label = paste0(gsub("_", " ", study), "\n(", disease, ")")
  )

label_order <- alpha_df |>
  dplyr::distinct(study, study_label) |>
  dplyr::arrange(study) |>
  dplyr::pull(study_label)

alpha_df <- alpha_df |>
  dplyr::mutate(study_label = factor(study_label, levels = label_order))

long_df <- alpha_df |>
  tidyr::pivot_longer(c(Observed, Chao1, Shannon, Simpson),
                      names_to  = "metric",
                      values_to = "value") |>
  dplyr::mutate(metric = factor(metric,
                                levels = c("Observed", "Chao1", "Shannon", "Simpson")))

# ---- Per-study Wilcoxon p for strip annotation --------------------------------
tests_study <- long_df |>
  dplyr::group_by(study_label, metric) |>
  dplyr::summarise(
    p = tryCatch(wilcox.test(value ~ case)$p.value,
                 error = function(e) NA_real_),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    label = dplyr::case_when(
      is.na(p)  ~ "",
      p < 0.001 ~ sprintf("p=%.2e", p),
      p < 0.05  ~ sprintf("p=%.3f", p),
      TRUE      ~ sprintf("p=%.2f",  p)
    )
  )

# ---- Plot --------------------------------------------------------------------
p <- ggplot(long_df, aes(x = case, y = value, fill = case)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.1, alpha = 0.7, colour = "grey30") +
  geom_text(data  = tests_study,
            aes(x = 1.5, y = Inf, label = label),
            vjust = 1.4, size = 2.8, inherit.aes = FALSE) +
  facet_grid(metric ~ study_label, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#AAAAAA", Case = "#CC6644")) +
  theme_bw(base_size = 10) +
  theme(
    strip.text.x    = element_text(size = 8),
    legend.position = "bottom",
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    panel.spacing.x = unit(0.3, "lines")
  ) +
  labs(
    x       = NULL,
    y       = "Alpha diversity value",
    fill    = "Group",
    caption = paste(
      "Studies ordered left-to-right: GERD (Ziganshina 2020, Qian 2023, Hao 2022, Kawar 2021)",
      "→ LPRD (Zheng 2024) → IBS-D (Li 2025, Tang 2023).",
      "Wilcoxon p shown per panel; significant (p < 0.05) values in bold."
    )
  ) +
  # Bold the significant p labels
  geom_text(data  = dplyr::filter(tests_study, !is.na(p), p < 0.05),
            aes(x = 1.5, y = Inf, label = label),
            vjust = 1.4, size = 2.8, fontface = "bold",
            inherit.aes = FALSE)

ggsave(file.path(FIGURES, "alpha_diversity.pdf"),
       p, width = 16, height = 8)
ggsave(file.path(FIGURES, "alpha_diversity.png"),
       p, width = 16, height = 8, dpi = 300)

message("Done. Saved:")
message("  figures/alpha_diversity.pdf  (replaces old unordered version)")
message("  figures/alpha_diversity.png")
