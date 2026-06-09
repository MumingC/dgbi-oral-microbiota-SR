# 04_differential_abundance.R — Per-study genus-level DA via MaAsLin2.
# Input: genus-collapsed SILVA counts. Outputs one MaAsLin2 folder per study.

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))

ps_list <- readRDS(file.path(CACHE, "ps_list.rds"))

for (s in studies$study) for (db in c("silva", "ehomd")) {
  ps <- ps_list[[paste0(s, "__", db)]]
  if (is.null(ps)) next
  message("MaAsLin2: ", s, " [", db, "]")

  ps_g <- phyloseq::tax_glom(ps, taxrank = "Genus", NArm = FALSE)

  # OTU matrix: rows = samples, cols = genera (MaAsLin2 convention).
  # Replace ASV IDs with genus names for readable output.
  tax <- as.data.frame(phyloseq::tax_table(ps_g))
  genus <- ifelse(is.na(tax$Genus) | tax$Genus == "",
                  paste0("Unclassified_", tax$Family), tax$Genus)
  otu_mat <- as(phyloseq::otu_table(ps_g), "matrix")
  if (phyloseq::taxa_are_rows(ps_g)) otu_mat <- t(otu_mat)  # -> samples × taxa
  colnames(otu_mat) <- make.unique(genus)
  otu_df <- as.data.frame(otu_mat, stringsAsFactors = FALSE)

  # Metadata: explicit S4 → data.frame, coerce group to character
  meta_raw <- as(phyloseq::sample_data(ps_g), "data.frame")
  grp_col <- grep("^group$", colnames(meta_raw), ignore.case = TRUE, value = TRUE)[1]
  meta_df <- data.frame(
    group = as.character(meta_raw[[grp_col]]),
    row.names = phyloseq::sample_names(ps_g),
    stringsAsFactors = FALSE
  )

  out_dir <- file.path(RESULTS, "maaslin2", paste0(s, "__", db))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  Maaslin2::Maaslin2(
    input_data     = otu_df,
    input_metadata = meta_df,
    output         = out_dir,
    fixed_effects  = "group",
    normalization  = "TSS",
    transform      = "LOG",
    min_prevalence = 0.1,
    min_abundance  = 0.0001,   # drop features present but vanishingly rare
    plot_heatmap   = TRUE,
    plot_scatter   = FALSE
  )
}

message("MaAsLin2 done. See results/maaslin2/<study>__{silva,ehomd}/all_results.tsv")
