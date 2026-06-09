# 02_alpha_diversity.R — Per-study alpha diversity + Wilcoxon tests.
# Shannon, observed features, Chao1 — computed on ASV-level SILVA tables.

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))

ps_list <- readRDS(file.path(CACHE, "ps_list.rds"))

alpha_df <- purrr::map_dfr(studies$study, function(s) {
  ps <- ps_list[[paste0(s, "__silva")]]
  if (is.null(ps)) return(NULL)
  er <- phyloseq::estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
  er$sample_id <- rownames(er)
  er$study <- s
  er$disease <- studies$disease[studies$study == s]
  er$group <- phyloseq::sample_data(ps)$group[match(er$sample_id, phyloseq::sample_names(ps))]
  er
})

# Keep GERD / LPRD / IBS-D separate (Zheng is driving the pooled GERD/LPRD signal alone)
alpha_df$disease_pool <- factor(alpha_df$disease, levels = c("GERD", "LPRD", "IBS-D"))
# Normalize group labels: control vs case (anything non-control)
alpha_df$case <- ifelse(tolower(alpha_df$group) %in% c("control", "healthy", "hc"),
                        "Control", "Case")

readr::write_tsv(alpha_df, file.path(RESULTS, "alpha_diversity.tsv"))

long_df <- alpha_df |>
  tidyr::pivot_longer(c(Observed, Chao1, Shannon, Simpson),
                      names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = factor(metric, levels = c("Observed", "Chao1", "Shannon", "Simpson")))

# Wilcoxon pooled per disease × metric (Case vs Control)
tests_pooled <- long_df |>
  dplyr::group_by(disease_pool, metric) |>
  dplyr::summarise(
    n_case    = sum(case == "Case"),
    n_control = sum(case == "Control"),
    p = tryCatch(wilcox.test(value ~ case)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  )
readr::write_tsv(tests_pooled, file.path(RESULTS, "alpha_wilcoxon_pooled.tsv"))

# Also keep per-study tests for sensitivity
tests_study <- long_df |>
  dplyr::group_by(study, metric) |>
  dplyr::summarise(
    p = tryCatch(wilcox.test(value ~ case)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  )
readr::write_tsv(tests_study, file.path(RESULTS, "alpha_wilcoxon_per_study.tsv"))

library(ggplot2)

# NOTE: the pooled-by-disease boxplot (alpha_diversity_pooled.{png,pdf}) was
# RETIRED on 2026-05-24. It pooled multiple studies' raw alpha values onto one
# axis with a per-disease Wilcoxon p that is confounded by between-study
# technical variation (e.g. IBS-D Observed p=5e-4, which is NOT significant
# under proper random-effects meta-analysis: pooled g=-0.73, 95% CI crossing 0,
# I2=85%). Inference now comes from the forest plot below. The per-study
# descriptive boxplot lives in 02b_alpha_diversity_by_disease.R.
# alpha_wilcoxon_pooled.tsv is still written above for the record only -- DO NOT
# cite its p-values in the manuscript.

# -------- Forest plot: per-study Hedges' g with random-effects pooling --------
# Effect size per study × metric: standardized mean difference (Case - Control).
# Positive = higher in Case (disease).

compute_smd <- function(df) {
  if (length(unique(df$case)) < 2) return(NULL)
  x_case <- df$value[df$case == "Case"]
  x_ctrl <- df$value[df$case == "Control"]
  n1 <- length(x_case); n2 <- length(x_ctrl)
  if (n1 < 2 || n2 < 2) return(NULL)
  es <- metafor::escalc(measure = "SMD",
                        m1i = mean(x_case), sd1i = sd(x_case), n1i = n1,
                        m2i = mean(x_ctrl), sd2i = sd(x_ctrl), n2i = n2)
  tibble::tibble(yi = es$yi, vi = es$vi, n_case = n1, n_ctrl = n2)
}

es_df <- long_df |>
  dplyr::group_by(study, disease_pool, metric) |>
  dplyr::group_modify(~ compute_smd(.x)) |>
  dplyr::ungroup()

# Random-effects meta-analysis per disease × metric
meta_rows <- es_df |>
  dplyr::group_by(disease_pool, metric) |>
  dplyr::group_modify(~{
    if (nrow(.x) < 2) return(tibble::tibble())
    fit <- tryCatch(metafor::rma(yi = .x$yi, vi = .x$vi, method = "DL"),
                    error = function(e) NULL)
    if (is.null(fit)) return(tibble::tibble())
    tibble::tibble(study = "Pooled (RE)",
                   yi = as.numeric(fit$beta), vi = fit$se^2,
                   ci.lb = fit$ci.lb, ci.ub = fit$ci.ub,
                   p = fit$pval, I2 = fit$I2, tau2 = fit$tau2,
                   k = nrow(.x), is_pooled = TRUE)
  }) |>
  dplyr::ungroup()

# Publication-standard forest layout:
#   * per-study Hedges' g as squares, area proportional to inverse-variance weight
#   * 95% CI as horizontal lines
#   * pooled DL random-effects estimate as a diamond whose horizontal width = 95% CI
#   * I^2 / k heterogeneity annotation beside each pooled diamond
# Note: escalc(measure = "SMD") above already applies the Hedges' g small-sample
# bias correction, so yi are Hedges' g (no extra correction needed).

metric_levels  <- c("Observed", "Chao1", "Shannon", "Simpson")
disease_levels <- c("GERD", "LPRD", "IBS-D")

study_rows <- es_df |>
  dplyr::mutate(
    ci.lb        = yi - 1.96 * sqrt(vi),
    ci.ub        = yi + 1.96 * sqrt(vi),
    is_pooled    = FALSE,
    weight       = 1 / vi,
    study        = as.character(study),
    disease_pool = as.character(disease_pool)
  ) |>
  dplyr::select(study, disease_pool, metric, yi, vi, ci.lb, ci.ub, is_pooled, weight)

pooled_rows <- meta_rows |>
  dplyr::mutate(
    study        = as.character(study),
    disease_pool = as.character(disease_pool),
    weight       = NA_real_
  ) |>
  dplyr::select(study, disease_pool, metric, yi, vi, ci.lb, ci.ub,
                is_pooled, weight, I2, tau2, k)

# Row order (top -> bottom): studies then pooled, within GERD -> LPRD -> IBS-D
layout <- dplyr::bind_rows(
  dplyr::mutate(dplyr::distinct(study_rows, disease_pool, study), is_pooled = FALSE),
  dplyr::mutate(dplyr::distinct(pooled_rows, disease_pool, study), is_pooled = TRUE)
) |>
  dplyr::mutate(disease_pool = factor(disease_pool, levels = disease_levels)) |>
  dplyr::arrange(disease_pool, is_pooled, study) |>
  dplyr::mutate(y = dplyr::n() - dplyr::row_number() + 1,
                disease_pool = as.character(disease_pool))

study_plot  <- dplyr::left_join(study_rows,  layout[, c("disease_pool", "study", "y")],
                                by = c("disease_pool", "study"))
pooled_plot <- dplyr::left_join(pooled_rows, layout[, c("disease_pool", "study", "y")],
                                by = c("disease_pool", "study"))

# Persist the underlying numbers (study SMDs + pooled estimates + heterogeneity)
readr::write_tsv(
  dplyr::bind_rows(
    dplyr::mutate(study_rows, I2 = NA_real_, tau2 = NA_real_, k = NA_integer_),
    pooled_rows
  ),
  file.path(RESULTS, "alpha_forest.tsv")
)

# Diamond polygons for the pooled rows (width = 95% CI, half-height = 0.25)
pooled_plot$gid <- paste(pooled_plot$disease_pool, pooled_plot$metric, sep = "__")
make_diamond <- function(df) {
  data.frame(
    disease_pool = df$disease_pool[1],
    metric       = df$metric[1],
    gid          = df$gid[1],
    x            = c(df$ci.lb[1], df$yi[1], df$ci.ub[1], df$yi[1]),
    y            = c(df$y[1], df$y[1] + 0.25, df$y[1], df$y[1] - 0.25)
  )
}
diamond_df <- purrr::map_dfr(split(pooled_plot, pooled_plot$gid), make_diamond)

het_df <- dplyr::mutate(pooled_plot,
                        label = sprintf("I²=%.0f%%, k=%d", I2, k))

# y-axis tick labels (studies + pooled); study underscores -> spaces
lab_df <- dplyr::mutate(layout,
                        label = ifelse(study == "Pooled (RE)",
                                       "Pooled (RE)", gsub("_", " ", study)))

# Factor ordering for facet rows (disease) and columns (metric)
fct <- function(d) dplyr::mutate(d,
  disease_pool = factor(disease_pool, levels = disease_levels),
  metric       = factor(metric, levels = metric_levels))
study_plot  <- fct(study_plot)
pooled_plot <- fct(pooled_plot)
diamond_df  <- fct(diamond_df)
het_df      <- fct(het_df)

forest_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(data = study_plot,
                aes(y = y, xmin = ci.lb, xmax = ci.ub),
                orientation = "y", width = 0, color = "grey30") +
  geom_point(data = study_plot,
             aes(x = yi, y = y, size = weight),
             shape = 15, color = "black") +
  geom_polygon(data = diamond_df,
               aes(x = x, y = y, group = gid),
               fill = "firebrick", color = "firebrick") +
  geom_text(data = het_df,
            aes(x = Inf, y = y, label = label),
            hjust = 1, vjust = -0.7, size = 2.5, color = "firebrick") +
  scale_size_area(max_size = 5, guide = "none") +
  scale_y_continuous(breaks = lab_df$y, labels = lab_df$label,
                     expand = ggplot2::expansion(add = 0.8)) +
  facet_grid(disease_pool ~ metric, scales = "free_y", space = "free_y") +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank()) +
  labs(x = "Standardized mean difference (Hedges' g, Case − Control)",
       y = NULL)
  # Caption removed 2026-05-26: figure legend now lives in the supplementary
  # text document (see Supplementary Figure S1 legend), per journal convention.

ggsave(file.path(FIGURES, "alpha_forest.png"),
       forest_plot, width = 12, height = 8, dpi = 300)
ggsave(file.path(FIGURES, "alpha_forest.pdf"),
       forest_plot, width = 12, height = 8)
print(forest_plot)

message("Alpha diversity done.")
message("  results/alpha_diversity.tsv — per-sample values")
message("  results/alpha_wilcoxon_pooled.tsv — Case vs Control per disease")
message("  results/alpha_wilcoxon_per_study.tsv — sensitivity per study")
message("  results/alpha_forest.tsv — per-study SMDs + REML pooled estimates")
message("  figures/alpha_forest.{png,pdf}")
