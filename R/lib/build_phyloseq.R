# lib/build_phyloseq.R — Construct phyloseq object from QIIME2 TSV exports.
#
# Expects per-study layout:
#   {study_dir}/metadata.tsv
#   {study_dir}/rooted-tree.nwk
#   {study_dir}/exported-{db}/asv-table.tsv
#   {study_dir}/exported-{db}/taxonomy.tsv
#
# db is "silva" or "ehomd".

build_phyloseq <- function(study_dir, db = c("silva", "ehomd")) {
  db <- match.arg(db)
  exp_dir <- file.path(study_dir, paste0("exported-", db))

  # ASV table: QIIME2 biom→tsv has "# Constructed from biom file" as line 1
  asv <- readr::read_tsv(file.path(exp_dir, "asv-table.tsv"),
                         skip = 1, show_col_types = FALSE)
  asv_mat <- as.matrix(asv[, -1])
  rownames(asv_mat) <- asv[[1]]
  storage.mode(asv_mat) <- "integer"

  # Taxonomy: QIIME2 outputs Feature ID / Taxon / Confidence
  tax <- readr::read_tsv(file.path(exp_dir, "taxonomy.tsv"),
                         show_col_types = FALSE)
  tax_split <- stringr::str_split_fixed(tax$Taxon, ";\\s*", 7)
  colnames(tax_split) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  tax_split[tax_split == ""] <- NA_character_
  rownames(tax_split) <- tax$`Feature ID`

  # Metadata
  meta <- readr::read_tsv(file.path(study_dir, "metadata.tsv"),
                          show_col_types = FALSE)
  meta_df <- as.data.frame(meta)
  rownames(meta_df) <- meta_df[[1]]

  # Tree (may not exist yet for Zheng/Tang post-merge — handle gracefully).
  # Binarize polytomies (phyloseq::UniFrac warns on multifurcating nodes; multi2di
  # adds zero-length branches to convert them to bifurcations — distances unchanged).
  tree_path <- file.path(study_dir, "rooted-tree.nwk")
  tree <- if (file.exists(tree_path)) {
    t <- ape::read.tree(tree_path)
    if (!ape::is.binary(t)) t <- ape::multi2di(t)
    t
  } else NULL

  # Align sample IDs across table + metadata
  common <- intersect(colnames(asv_mat), rownames(meta_df))
  if (!length(common)) stop("No sample IDs in common between ASV table and metadata for ", study_dir)
  asv_mat <- asv_mat[, common, drop = FALSE]
  meta_df <- meta_df[common, , drop = FALSE]

  args <- list(
    phyloseq::otu_table(asv_mat, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax_split),
    phyloseq::sample_data(meta_df)
  )
  if (!is.null(tree)) args <- c(args, list(phyloseq::phy_tree(tree)))
  do.call(phyloseq::phyloseq, args)
}
