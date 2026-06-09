# qc_silva_ehomd_crossref.R — Cross-reference genus-level taxonomy between
# SILVA v138 and eHOMD v16.02 for all 7 studies.
#
# Outputs:
#   results/taxonomy_crossref.tsv        — all genus-level disagreements
#   results/taxonomy_crossref_summary.tsv — aggregated genus mapping table
#   Console messages for top disagreements and the specific
#   Pseudopropionibacterium (SILVA) vs Arachnia (eHOMD) check.

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))

# Load cached phyloseq list
ps_list <- readRDS(file.path(CACHE, "ps_list.rds"))

# Strip QIIME2 rank prefix (e.g. "g__Streptococcus" -> "Streptococcus")
strip_prefix <- function(x) {
  out <- stringr::str_remove(x, "^[dkpcofgs]__")
  out[is.na(out) | out == "" | out == "__"] <- NA_character_
  out
}

# Loop through studies: extract and join genus assignments
all_disagreements <- list()

for (s in studies$study) {
  key_silva <- paste0(s, "__silva")
  key_ehomd <- paste0(s, "__ehomd")

  ps_s <- ps_list[[key_silva]]
  ps_e <- ps_list[[key_ehomd]]

  if (is.null(ps_s) || is.null(ps_e)) {
    message("Skipping ", s, " — missing phyloseq object")
    next
  }

  # Extract tax tables as data frames, keep Genus column
  tax_s <- as.data.frame(phyloseq::tax_table(ps_s), stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("ASV") %>%
    transmute(ASV,
              Genus_silva = strip_prefix(Genus),
              Family_silva = strip_prefix(Family),
              Phylum_silva = strip_prefix(Phylum))

  tax_e <- as.data.frame(phyloseq::tax_table(ps_e), stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("ASV") %>%
    transmute(ASV,
              Genus_ehomd = strip_prefix(Genus),
              Family_ehomd = strip_prefix(Family),
              Phylum_ehomd = strip_prefix(Phylum))

  # Inner join on ASV ID (should be the same set for both databases)
  joined <- inner_join(tax_s, tax_e, by = "ASV") %>%
    mutate(study = s)

  # Flag disagreements (both assigned at genus level but different names)
  joined <- joined %>%
    mutate(
      both_assigned = !is.na(Genus_silva) & !is.na(Genus_ehomd),
      agree         = both_assigned & (Genus_silva == Genus_ehomd),
      disagree      = both_assigned & (Genus_silva != Genus_ehomd),
      silva_only    = !is.na(Genus_silva) & is.na(Genus_ehomd),
      ehomd_only    = is.na(Genus_silva) & !is.na(Genus_ehomd),
      both_na       = is.na(Genus_silva) & is.na(Genus_ehomd)
    )

  # Summary for this study
  n_total   <- nrow(joined)
  n_agree   <- sum(joined$agree, na.rm = TRUE)
  n_disagr  <- sum(joined$disagree, na.rm = TRUE)
  n_silva   <- sum(joined$silva_only, na.rm = TRUE)
  n_ehomd   <- sum(joined$ehomd_only, na.rm = TRUE)
  n_bothna  <- sum(joined$both_na, na.rm = TRUE)

  message(sprintf(
    "%s: %d ASVs | agree=%d | disagree=%d | SILVA-only=%d | eHOMD-only=%d | both-NA=%d",
    s, n_total, n_agree, n_disagr, n_silva, n_ehomd, n_bothna
  ))

  all_disagreements[[s]] <- joined %>% filter(disagree)
}

# Combine all disagreements across studies
disagreements <- bind_rows(all_disagreements)

message("\n=== Total genus-level disagreements across all studies: ",
        nrow(disagreements), " ASVs ===\n")

# Write full disagreement table
write_tsv(
  disagreements %>%
    select(study, ASV, Genus_silva, Genus_ehomd,
           Family_silva, Family_ehomd, Phylum_silva, Phylum_ehomd),
  file.path(RESULTS, "taxonomy_crossref.tsv")
)
message("Wrote: ", file.path(RESULTS, "taxonomy_crossref.tsv"))

# Aggregated genus mapping table (SILVA -> eHOMD)
genus_map <- disagreements %>%
  count(Genus_silva, Genus_ehomd, name = "n_asvs") %>%
  arrange(desc(n_asvs))

write_tsv(genus_map, file.path(RESULTS, "taxonomy_crossref_summary.tsv"))
message("Wrote: ", file.path(RESULTS, "taxonomy_crossref_summary.tsv"))

# Print top 30 disagreement pairs
message("\n--- Top 30 genus-level disagreements (SILVA -> eHOMD) ---")
print(as.data.frame(head(genus_map, 30)), row.names = FALSE)

# Which studies contributed to each mapping
genus_map_studies <- disagreements %>%
  group_by(Genus_silva, Genus_ehomd) %>%
  summarise(
    n_asvs  = n(),
    studies = paste(unique(study), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_asvs))

message("\n--- Genus mapping with study sources ---")
print(as.data.frame(head(genus_map_studies, 30)), row.names = FALSE)

# Specific check: Pseudopropionibacterium (SILVA) vs Arachnia (eHOMD)
message("\n=== Specific check: Pseudopropionibacterium (SILVA) vs Arachnia (eHOMD) ===")

pseudo_arachnia <- disagreements %>%
  filter(Genus_silva == "Pseudopropionibacterium" & Genus_ehomd == "Arachnia")

if (nrow(pseudo_arachnia) > 0) {
  message("FOUND: ", nrow(pseudo_arachnia),
          " ASVs classified as Pseudopropionibacterium by SILVA and Arachnia by eHOMD")
  message("Studies: ", paste(unique(pseudo_arachnia$study), collapse = ", "))
  print(as.data.frame(pseudo_arachnia %>%
    select(study, ASV, Genus_silva, Genus_ehomd, Family_silva, Family_ehomd)))
} else {
  message("NOT FOUND among disagreements.")
}

# Also check where Pseudopropionibacterium appears in SILVA regardless of eHOMD
pseudo_all <- bind_rows(all_disagreements) # already only disagreements
# Rebuild from full joined data to include agreements too
pseudo_full <- list()
for (s in studies$study) {
  ps_s <- ps_list[[paste0(s, "__silva")]]
  ps_e <- ps_list[[paste0(s, "__ehomd")]]
  if (is.null(ps_s) || is.null(ps_e)) next

  tax_s <- as.data.frame(phyloseq::tax_table(ps_s), stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("ASV") %>%
    transmute(ASV, Genus_silva = strip_prefix(Genus))

  tax_e <- as.data.frame(phyloseq::tax_table(ps_e), stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("ASV") %>%
    transmute(ASV, Genus_ehomd = strip_prefix(Genus))

  joined <- inner_join(tax_s, tax_e, by = "ASV") %>%
    mutate(study = s)

  hits <- joined %>% filter(Genus_silva == "Pseudopropionibacterium" |
                            Genus_ehomd == "Arachnia")
  if (nrow(hits) > 0) pseudo_full[[s]] <- hits
}
pseudo_full_df <- bind_rows(pseudo_full)

if (nrow(pseudo_full_df) > 0) {
  message("\n--- All ASVs with Pseudopropionibacterium (SILVA) or Arachnia (eHOMD) ---")
  print(as.data.frame(pseudo_full_df))
} else {
  message("No ASVs found with either Pseudopropionibacterium (SILVA) or Arachnia (eHOMD) in any study.")
}

# Per-study agreement rate summary
message("\n=== Per-study agreement rate (genus level, both databases assigned) ===")
rate_summary <- list()
for (s in studies$study) {
  ps_s <- ps_list[[paste0(s, "__silva")]]
  ps_e <- ps_list[[paste0(s, "__ehomd")]]
  if (is.null(ps_s) || is.null(ps_e)) next

  tax_s <- as.data.frame(phyloseq::tax_table(ps_s), stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("ASV") %>%
    transmute(ASV, Genus_silva = strip_prefix(Genus))

  tax_e <- as.data.frame(phyloseq::tax_table(ps_e), stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("ASV") %>%
    transmute(ASV, Genus_ehomd = strip_prefix(Genus))

  joined <- inner_join(tax_s, tax_e, by = "ASV")
  both   <- joined %>% filter(!is.na(Genus_silva) & !is.na(Genus_ehomd))
  n_both <- nrow(both)
  n_agree <- sum(both$Genus_silva == both$Genus_ehomd)

  rate_summary[[s]] <- tibble(
    study        = s,
    total_asvs   = nrow(joined),
    both_assigned = n_both,
    agree        = n_agree,
    disagree     = n_both - n_agree,
    pct_agree    = round(100 * n_agree / n_both, 1)
  )
}
rate_df <- bind_rows(rate_summary)
print(as.data.frame(rate_df), row.names = FALSE)

write_tsv(rate_df, file.path(RESULTS, "taxonomy_crossref_rates.tsv"))
message("\nWrote: ", file.path(RESULTS, "taxonomy_crossref_rates.tsv"))
message("Done.")
