# 03_beta_diversity.R — Per-study beta diversity on compositional-transformed data.
# Bray-Curtis + weighted/unweighted UniFrac → PCoA, PERMANOVA (adonis2), PERMADISP (betadisper).
# Convention matches supervisor's Qiime2.Rmd: microbiome::transform("compositional") before distance.

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))

library(ggplot2)
library(vegan)

set.seed(123)  # reproducible adonis2 / betadisper permutation p-values

ps_list <- readRDS(file.path(CACHE, "ps_list.rds"))

permanova_rows <- list()
permdisp_rows  <- list()
pcoa_rows      <- list()

for (s in studies$study) {
  ps <- ps_list[[paste0(s, "__silva")]]
  if (is.null(ps)) next
  message("Beta: ", s)

  # Compositional transform (relative abundance)
  ps_rel <- microbiome::transform(ps, "compositional")
  meta   <- data.frame(phyloseq::sample_data(ps_rel))
  meta$case <- ifelse(tolower(meta$group) %in% c("control", "healthy", "hc"),
                      "Control", "Case")

  for (metric in c("bray", "wunifrac", "unifrac")) {
    d <- tryCatch(phyloseq::distance(ps_rel, method = metric),
                  error = function(e) NULL)
    if (is.null(d)) next

    # PERMANOVA
    ad <- tryCatch(vegan::adonis2(d ~ case, data = meta, permutations = 999),
                   error = function(e) NULL)
    if (!is.null(ad)) {
      permanova_rows[[length(permanova_rows) + 1]] <- tibble::tibble(
        study = s, disease = studies$disease[studies$study == s],
        metric = metric, R2 = ad$R2[1], F = ad$F[1], p = ad$`Pr(>F)`[1]
      )
    }

    # PERMADISP — dispersion homogeneity check
    bd <- tryCatch(vegan::betadisper(d, meta$case), error = function(e) NULL)
    if (!is.null(bd)) {
      pd <- vegan::permutest(bd, permutations = 999)
      permdisp_rows[[length(permdisp_rows) + 1]] <- tibble::tibble(
        study = s, disease = studies$disease[studies$study == s],
        metric = metric, F = pd$tab$F[1], p = pd$tab$`Pr(>F)`[1]
      )
    }

    # PCoA coordinates for plotting (Bray only — most commonly reported)
    if (metric == "bray") {
      ord <- phyloseq::ordinate(ps_rel, method = "PCoA", distance = d)
      axes <- as.data.frame(ord$vectors[, 1:2])
      colnames(axes) <- c("Axis1", "Axis2")
      eigs <- ord$values$Relative_eig[1:2]
      axes$sample_id <- rownames(axes)
      axes$study <- s
      axes$disease <- studies$disease[studies$study == s]
      axes$case <- meta$case[match(axes$sample_id, rownames(meta))]
      axes$var_axis1 <- round(eigs[1] * 100, 1)
      axes$var_axis2 <- round(eigs[2] * 100, 1)
      pcoa_rows[[length(pcoa_rows) + 1]] <- axes
    }
  }
}

permanova_df <- dplyr::bind_rows(permanova_rows)
permdisp_df  <- dplyr::bind_rows(permdisp_rows)
pcoa_df      <- dplyr::bind_rows(pcoa_rows)

readr::write_tsv(permanova_df, file.path(RESULTS, "permanova.tsv"))
readr::write_tsv(permdisp_df,  file.path(RESULTS, "permdisp.tsv"))
readr::write_tsv(pcoa_df,      file.path(RESULTS, "pcoa_bray.tsv"))

# -------- PCoA plot: single facet grid, GERD+LPRD row on top, IBS-D row below --------
pcoa_annot <- permanova_df |>
  dplyr::filter(metric == "bray") |>
  dplyr::mutate(label = sprintf("R²=%.3f, p=%.3f", R2, p))

# Study order: GERD (4) → LPRD (1) → IBS-D (2), so row 1 is 5 panels, row 2 is 2.
study_order <- c("Ziganshina_2020", "Qian_2023", "Hao_2022", "Kawar_2021",
                 "Zheng_2024", "Li_2025", "Tang_2023")
study_order <- intersect(study_order, unique(pcoa_df$study))

pcoa_df$study_lab <- factor(
  paste0(pcoa_df$study, "\n(", pcoa_df$disease, ")"),
  levels = paste0(study_order,
                  "\n(", studies$disease[match(study_order, studies$study)], ")")
)
pcoa_annot$study_lab <- factor(
  paste0(pcoa_annot$study, "\n(", pcoa_annot$disease, ")"),
  levels = levels(pcoa_df$study_lab)
)

pcoa_plot <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = case)) +
  geom_point(size = 1.5, alpha = 0.8) +
  stat_ellipse(aes(group = case), level = 0.8, linewidth = 0.5) +
  geom_text(data = pcoa_annot,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.3, inherit.aes = FALSE, size = 3) +
  facet_wrap(~ study_lab, scales = "free", ncol = 5) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9, lineheight = 0.9)) +
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2", color = "Group")
  # Caption removed 2026-05-26: figure legend now lives in the supplementary
  # legends document (Product\Supplementary_Figure_Legends.md, Figure S2).

ggsave(file.path(FIGURES, "pcoa_bray.png"), pcoa_plot, width = 15, height = 7, dpi = 300)
ggsave(file.path(FIGURES, "pcoa_bray.pdf"), pcoa_plot, width = 15, height = 7)
print(pcoa_plot)

# -------- Forest plot of PERMANOVA R² --------
# R² is a per-study effect size; plot with 95% CI unavailable from adonis2 by default, so show point only
perm_forest <- permanova_df |>
  dplyr::mutate(
    disease = factor(disease, levels = c("GERD", "LPRD", "IBS-D")),
    sig = ifelse(p < 0.05, "p<0.05", "ns")
  )

permanova_plot <- ggplot(perm_forest, aes(x = R2, y = study, color = sig)) +
  geom_point(size = 3) +
  # Place p-value label just left of the point so it never falls off the panel
  geom_text(aes(label = sprintf("p=%.3f", p)),
            hjust = 1.15, size = 3, color = "black") +
  scale_color_manual(values = c("p<0.05" = "firebrick", "ns" = "grey50")) +
  scale_x_continuous(expand = expansion(mult = c(0.25, 0.1))) +
  facet_grid(disease ~ metric, scales = "free_y", space = "free_y") +
  theme_bw(base_size = 10) +
  labs(x = "PERMANOVA R² (Case vs Control)", y = NULL, color = NULL)
  # Caption removed 2026-05-26: figure legend now lives in the supplementary
  # legends document (Product\Supplementary_Figure_Legends.md, Figure S3).

ggsave(file.path(FIGURES, "permanova_forest.png"), permanova_plot, width = 11, height = 8, dpi = 300)
ggsave(file.path(FIGURES, "permanova_forest.pdf"), permanova_plot, width = 11, height = 8)
print(permanova_plot)

message("Beta diversity done.")
message("  results/permanova.tsv — Case vs Control R², p")
message("  results/permdisp.tsv  — dispersion test (validates PERMANOVA)")
message("  results/pcoa_bray.tsv — PCoA coordinates")
message("  figures/pcoa_bray.{png,pdf}")
message("  figures/permanova_forest.{png,pdf}")
