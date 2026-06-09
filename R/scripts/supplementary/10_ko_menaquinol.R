# 10_ko_menaquinol.R — Supplementary KO-level heatmap for menaquinol biosynthesis.
# Companion to 09_heatmap_picrust2.R: verifies the MetaCyc menaquinol pathway
# signal (PWY-5845/5850/5860/5862/5896, all padj<0.05 ↓ in IBS-D) at the
# gene-family (KEGG KO) resolution.
#
# Narrative: at KO-level, only K19222 (MenI, DHNA-CoA hydrolase) reaches
# meta padj<0.05 in IBS-D. The classical men genes (menA-H) and the
# alternative futalosine mqn genes show weak-but-concordant reductions that
# individually don't cross the gene-level multiple-testing threshold. This is
# how MinPath aggregates weak coherent gene signals into a significant pathway
# result — and this figure makes that coherence visible.
#
# Row filter (hard-coded): curated menaquinone KO set in
# lib/menaquinol_ko_names.tsv. No padj/concordance filter — we show all genes
# in the pathway so the reader can see the diffuse pattern.
#
# Outputs (figures/heatmap/):
#   heatmap_ko_menaquinol_by_padj.pdf   # rows ordered by min meta padj
#   heatmap_ko_menaquinol_by_clust.pdf  # rows ordered by hclust on coef

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))
source(here::here("lib", "heatmap_helpers.R"))

install_heatmap_deps()

thresholds <- list(star1 = 0.25, star2 = 0.05, meta_sig = 0.05, color_limit = 3)

HEATDIR <- file.path(FIGURES, "heatmap")
dir.create(HEATDIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load the curated KO->gene name map ------------------------------------
ko_map_path <- here::here("lib", "menaquinol_ko_names.tsv")
if (!file.exists(ko_map_path)) {
  stop("KO name map not found at ", ko_map_path)
}
ko_map <- readr::read_tsv(ko_map_path, show_col_types = FALSE)
menaquinol_kos <- ko_map$ko_id
message("Menaquinol KO set: ", length(menaquinol_kos), " KOs (",
        sum(ko_map$pathway_branch == "classical"), " classical + ",
        sum(ko_map$pathway_branch == "futalosine"), " futalosine)")

# Row filter: keep only KOs in the curated menaquinol set that are also
# present in the meta table (PICRUSt2 may not have tested a KO if it had
# too many zeros).
row_filter_menaquinol <- function(meta_df, per_study_df = NULL) {
  present <- intersect(menaquinol_kos, unique(meta_df$feature))
  missing <- setdiff(menaquinol_kos, present)
  if (length(missing)) {
    message("  note: ", length(missing),
            " menaquinol KO(s) not in meta table (dropped): ",
            paste(missing, collapse = ", "))
  }
  present
}

# Row labels: "K19222 — menI" style, with pathway-branch suffix so the reader
# can spot classical vs futalosine.
ko_label_fn <- function(ids) {
  lookup <- setNames(
    sprintf("%s \u2014 %s (%s)",
            ko_map$ko_id, ko_map$gene,
            ifelse(ko_map$pathway_branch == "classical", "classical", "futalosine")),
    ko_map$ko_id
  )
  hit <- lookup[ids]
  ifelse(is.na(hit), ids, unname(hit))
}

message("\n=== PICRUSt2 KO \u2014 menaquinol biosynthesis ===")

dat <- build_heatmap_data(
  per_study_path = file.path(RESULTS, "picrust2_all_studies__ko.tsv"),
  meta_path      = file.path(RESULTS, "meta_analysis_picrust2__ko.tsv"),
  row_filter_fn  = row_filter_menaquinol
)
n_rows <- nrow(dat$coef)
message("  rows: ", n_rows)
if (n_rows == 0) stop("Nothing to plot.")

# Small figure — 14 KOs -> modest height.
height <- min(30, max(5, 0.32 * n_rows + 3))
for (order_mode in c("padj", "clust")) {
  ht <- build_heatmap(
    dat,
    order_mode   = order_mode,
    title_prefix = sprintf(
      "PICRUSt2 KO \u2014 menaquinol biosynthesis (curated, %d KOs)",
      n_rows
    ),
    row_label_fn = ko_label_fn,
    thresholds   = thresholds
  )
  out <- file.path(HEATDIR,
                   sprintf("heatmap_ko_menaquinol_by_%s.pdf", order_mode))
  pdf(out, width = 14, height = height)
  ComplexHeatmap::draw(
    ht,
    merge_legend           = TRUE,
    heatmap_legend_side    = "right",
    annotation_legend_side = "right"
  )
  dev.off()
  message("  wrote ", out)
}

message("\nDone. Menaquinol KO heatmap in ", HEATDIR)
