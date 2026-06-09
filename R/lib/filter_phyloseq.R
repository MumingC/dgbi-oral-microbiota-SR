# lib/filter_phyloseq.R — Standard QC filter (matches supervisor's Qiime2.Rmd).
# Remove Eukaryota / Chloroplast / unknown Kingdom, prune low-read samples, drop empty taxa.

filter_phyloseq <- function(ps, min_reads = 2000) {
  ps <- phyloseq::subset_taxa(ps, !Kingdom %in% c("Eukaryota", NA))
  ps <- phyloseq::subset_taxa(ps, is.na(Order) | Order != "Chloroplast")
  ps <- phyloseq::prune_samples(phyloseq::sample_sums(ps) >= min_reads, ps)
  ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
  ps
}
