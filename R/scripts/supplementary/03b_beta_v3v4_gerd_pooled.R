# 03b_beta_v3v4_gerd_pooled.R
# Pooled beta diversity across the 3 V3-V4 GERD studies: Ziganshina, Qian, Hao.
# Approach: collapse each study to genus → merge genus tables → Bray-Curtis → PCoA.
# Handles the "ASV IDs differ per study" problem by using genus names as the common row key.

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))

library(ggplot2)
library(vegan)

set.seed(123)  # reproducible adonis2 / betadisper permutation p-values

ps_list <- readRDS(file.path(CACHE, "ps_list.rds"))

target_studies <- c("Ziganshina_2020", "Qian_2023", "Hao_2022")  # V3-V4 GERD

# ---- Build per-study genus-level tables (rows = Genus, cols = samples) ----
genus_tables <- list()
meta_list    <- list()
for (s in target_studies) {
  ps <- ps_list[[paste0(s, "__silva")]]
  if (is.null(ps)) { message("Skip (missing): ", s); next }
  ps_g <- phyloseq::tax_glom(ps, taxrank = "Genus", NArm = FALSE)

  # Replace ASV IDs with genus names for the row key
  tax <- as.data.frame(phyloseq::tax_table(ps_g))
  genus <- ifelse(is.na(tax$Genus) | tax$Genus == "",
                  paste0("Unclassified_", tax$Family), tax$Genus)
  otu <- as(phyloseq::otu_table(ps_g), "matrix")
  if (!phyloseq::taxa_are_rows(ps_g)) otu <- t(otu)
  rownames(otu) <- make.unique(genus)
  genus_tables[[s]] <- otu

  # phyloseq sample_data is S4; convert explicitly to a real data.frame.
  meta_raw <- as(phyloseq::sample_data(ps_g), "data.frame")
  sid <- phyloseq::sample_names(ps_g)
  # Find the 'group' column case-insensitively
  grp_col <- grep("^group$", colnames(meta_raw), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(grp_col)) stop("No 'group' column found in metadata for ", s)
  grp_vals <- as.character(meta_raw[[grp_col]])

  meta <- data.frame(
    sample_id = paste0(s, "__", sid),
    study     = s,
    disease   = studies$disease[studies$study == s],
    group     = grp_vals,
    case      = ifelse(tolower(grp_vals) %in% c("control", "healthy", "hc"),
                       "Control", "Case"),
    row.names = paste0(s, "__", sid),
    stringsAsFactors = FALSE
  )
  colnames(genus_tables[[s]]) <- paste0(s, "__", colnames(genus_tables[[s]]))
  meta_list[[s]] <- meta
}

# ---- Merge: union of genera across studies, missing = 0 ----
all_genera <- unique(unlist(lapply(genus_tables, rownames)))
merged_mat <- matrix(0, nrow = length(all_genera),
                     ncol = sum(sapply(genus_tables, ncol)))
rownames(merged_mat) <- all_genera
colnames(merged_mat) <- unlist(lapply(genus_tables, colnames))
for (s in names(genus_tables)) {
  m <- genus_tables[[s]]
  merged_mat[rownames(m), colnames(m)] <- m
}
# unname() to stop rbind prepending list names (e.g. "Ziganshina_2020.<id>") to rownames
merged_meta <- do.call(rbind, unname(meta_list))
merged_meta <- merged_meta[colnames(merged_mat), , drop = FALSE]

message("Merged: ", nrow(merged_mat), " genera × ", ncol(merged_mat), " samples")
message("  Studies: ", paste(unique(merged_meta$study), collapse = ", "))
message("  Case: ",    sum(merged_meta$case == "Case", na.rm = TRUE),
        " | Control: ", sum(merged_meta$case == "Control", na.rm = TRUE))

# ---- Relative abundance (compositional) + Bray-Curtis ----
rel_mat <- apply(merged_mat, 2, function(x) if (sum(x) > 0) x / sum(x) else x)
d_bray <- vegan::vegdist(t(rel_mat), method = "bray")

# ---- PERMANOVA: naive + strata (within-study) ----
ad_naive <- vegan::adonis2(d_bray ~ case, data = merged_meta, permutations = 999)
ad_strat <- vegan::adonis2(d_bray ~ case, data = merged_meta,
                           strata = merged_meta$study, permutations = 999)
ad_study <- vegan::adonis2(d_bray ~ study + case, data = merged_meta, permutations = 999,
                           by = "margin")

permanova_pooled <- tibble::tibble(
  model = c("case (naive)", "case | study (strata)", "study + case (adjusted)"),
  R2_case = c(ad_naive$R2[1], ad_strat$R2[1], ad_study$R2[2]),
  F_case  = c(ad_naive$F[1],  ad_strat$F[1],  ad_study$F[2]),
  p_case  = c(ad_naive$`Pr(>F)`[1], ad_strat$`Pr(>F)`[1], ad_study$`Pr(>F)`[2])
)
readr::write_tsv(permanova_pooled, file.path(RESULTS, "permanova_v3v4_gerd.tsv"))

# ---- PERMADISP ----
bd <- vegan::betadisper(d_bray, merged_meta$case)
pd <- vegan::permutest(bd, permutations = 999)
message("PERMADISP p = ", signif(pd$tab$`Pr(>F)`[1], 3),
        " (if significant, PERMANOVA may reflect dispersion, not location shift)")

# ---- PCoA ----
pcoa <- ape::pcoa(d_bray)
eigs <- pcoa$values$Relative_eig[1:2]
coords <- data.frame(Axis1 = pcoa$vectors[, 1], Axis2 = pcoa$vectors[, 2])
coords$sample_id <- rownames(coords)
coords <- cbind(coords, merged_meta[coords$sample_id, c("study", "case")])
readr::write_tsv(coords, file.path(RESULTS, "pcoa_v3v4_gerd.tsv"))

lab_r2 <- sprintf("Naive: R²=%.3f, p=%.3f\nStrata by study: R²=%.3f, p=%.3f",
                  ad_naive$R2[1], ad_naive$`Pr(>F)`[1],
                  ad_strat$R2[1], ad_strat$`Pr(>F)`[1])

pcoa_plot <- ggplot(coords, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = case, shape = study), size = 2.5, alpha = 0.85) +
  stat_ellipse(aes(group = case, color = case), level = 0.8, linewidth = 0.5) +
  annotate("text", x = -Inf, y = Inf, label = lab_r2,
           hjust = -0.05, vjust = 1.2, size = 3.3) +
  theme_bw(base_size = 11) +
  labs(x = sprintf("PCoA Axis 1 (%.1f%%)", eigs[1] * 100),
       y = sprintf("PCoA Axis 2 (%.1f%%)", eigs[2] * 100),
       color = "Group", shape = "Study")
  # Title and caption removed 2026-05-26: figure legend now lives in the
  # supplementary legends document (Product\Supplementary_Figure_Legends.md,
  # Figure S2 after renumbering). The lab_r2 annotation inside the plot panel
  # is retained because it carries the actual PERMANOVA statistics.

ggsave(file.path(FIGURES, "pcoa_v3v4_gerd.png"), pcoa_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(FIGURES, "pcoa_v3v4_gerd.pdf"), pcoa_plot, width = 8, height = 6)
print(pcoa_plot)

message("\nPooled V3-V4 GERD beta diversity done.")
message("  results/permanova_v3v4_gerd.tsv")
message("  results/pcoa_v3v4_gerd.tsv")
message("  figures/pcoa_v3v4_gerd.{png,pdf}")
