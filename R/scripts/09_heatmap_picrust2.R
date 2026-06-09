# 09_heatmap_picrust2.R — Duvallet-style effect-size heatmap for PICRUSt2
# MetaCyc pathways. Companion to 08_heatmap.R; same layout/conventions,
# driven by the shared lib/heatmap_helpers.R.
#
# Row filter (default): pathways with meta padj < MAIN_FILTER_PADJ in at least
# one disease. At padj<0.05 this is 22 pathways, all IBS-D (GERD has 0 hits).
# This is the honest evidence of functional convergence in IBS-D that the
# taxa meta-analysis lacks.
#
# Outputs (figures/heatmap/):
#   heatmap_picrust2_pathway_by_padj.pdf   # rows ordered by min meta padj
#   heatmap_picrust2_pathway_by_clust.pdf  # rows ordered by hclust on coef
#
# To loosen the filter:
#   MAIN_FILTER_PADJ <- 0.25  -> 38 pathways incl. some GERD
#   MAIN_FILTER_PADJ <- 1     -> all 377 tested pathways

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))
source(here::here("lib", "heatmap_helpers.R"))

install_heatmap_deps()

# ---- Top-level knobs --------------------------------------------------------
MAIN_FILTER_PADJ <- 0.05
thresholds <- list(star1 = 0.25, star2 = 0.05, meta_sig = 0.05, color_limit = 3)

HEATDIR <- file.path(FIGURES, "heatmap")
dir.create(HEATDIR, recursive = TRUE, showWarnings = FALSE)

# Row filter: meta padj < MAIN_FILTER_PADJ in at least one disease.
# Signature takes (meta_df, per_study_df) for API symmetry with
# row_filter_concordance() in lib/heatmap_helpers.R; per_study_df is unused here.
row_filter_sig_padj <- function(meta_df, per_study_df = NULL) {
  meta_df %>%
    dplyr::filter(!is.na(padj), padj < MAIN_FILTER_PADJ) %>%
    dplyr::distinct(feature) %>%
    dplyr::pull(feature)
}

# Row labels: map MetaCyc codes -> short English names via lib/metacyc_pathway_names.tsv.
# Fallback to raw ID if a code is missing from the map (so broadening the filter
# doesn't silently drop labels).
load_metacyc_name_fn <- function(
  path = here::here("lib", "metacyc_pathway_names.tsv")
) {
  if (!file.exists(path)) {
    warning("MetaCyc name map not found at ", path, "; keeping raw IDs.")
    return(identity)
  }
  m <- readr::read_tsv(path, show_col_types = FALSE)
  nm <- setNames(m$description, m$pathway_id)
  function(ids) {
    hit <- nm[ids]
    ifelse(is.na(hit), ids, unname(hit))
  }
}
metacyc_name_fn <- load_metacyc_name_fn()

# ---- Right-side functional category annotation ------------------------------
# Maps MetaCyc pathway IDs -> category label for the right-side colour bar.
# Unlisted pathways fall through to "Other metabolism".
PATHWAY_CATEGORIES <- c(
  # Menaquinol / demethylmenaquinol biosynthesis (5 pathways)
  "PWY-5845"           = "Menaquinol biosynthesis",
  "PWY-5850"           = "Menaquinol biosynthesis",
  "PWY-5860"           = "Menaquinol biosynthesis",
  "PWY-5862"           = "Menaquinol biosynthesis",
  "PWY-5896"           = "Menaquinol biosynthesis",
  # Aromatic amino acid / phenylpropanoate metabolism (3 pathways)
  "ALL-CHORISMATE-PWY" = "Aromatic amino acid\nmetabolism",
  "TYRFUMCAT-PWY"      = "Aromatic amino acid\nmetabolism",
  "P281-PWY"           = "Aromatic amino acid\nmetabolism",
  # LPS / cell wall biosynthesis — enriched in IBS-D (4 pathways)
  "NAGLIPASYN-PWY"     = "LPS / cell wall\nbiosynthesis",
  "PWY-1269"           = "LPS / cell wall\nbiosynthesis",
  "PWY-6467"           = "LPS / cell wall\nbiosynthesis",
  "PWY-6471"           = "LPS / cell wall\nbiosynthesis"
  # All remaining pathways are labelled "Other metabolism" via the fallback below.
)

CATEGORY_COLORS <- c(
  "Menaquinol biosynthesis"        = "#e07b39",   # warm orange
  "Aromatic amino acid\nmetabolism" = "#7b5ea7",  # purple
  "LPS / cell wall\nbiosynthesis"  = "#c0392b",   # red (enriched)
  "Other metabolism"               = "grey82"
)

# Build a rowAnnotation aligned to the supplied pathway IDs (original matrix row order).
# ComplexHeatmap reorders annotation rows together with the main heatmap.
make_category_ann <- function(row_ids) {
  cats <- dplyr::case_when(
    row_ids %in% names(PATHWAY_CATEGORIES) ~ PATHWAY_CATEGORIES[row_ids],
    TRUE                                   ~ "Other metabolism"
  )
  ComplexHeatmap::rowAnnotation(
    Category = cats,
    col      = list(Category = CATEGORY_COLORS),
    width    = grid::unit(3, "mm"),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Category = list(
        title     = "Pathway category",
        labels_gp = grid::gpar(fontsize = 8),
        title_gp  = grid::gpar(fontsize = 9, fontface = "bold")
      )
    )
  )
}

# Report unmapped features so the name map can be extended if filter broadens.
check_unmapped <- function(ids, fn) {
  mapped <- fn(ids)
  unmapped <- ids[mapped == ids]
  if (length(unmapped)) {
    message("  note: ", length(unmapped),
            " pathway(s) have no English name in the map (showing raw ID): ",
            paste(head(unmapped, 6), collapse = ", "),
            if (length(unmapped) > 6) ", ..." else "")
  }
}

feature_type <- "pathway"
message("\n=== PICRUSt2 ", feature_type,
        " (padj < ", MAIN_FILTER_PADJ, ") ===")

dat <- build_heatmap_data(
  per_study_path = file.path(RESULTS,
                             paste0("picrust2_all_studies__", feature_type, ".tsv")),
  meta_path      = file.path(RESULTS,
                             paste0("meta_analysis_picrust2__", feature_type, ".tsv")),
  row_filter_fn  = row_filter_sig_padj
)
n_rows <- nrow(dat$coef)
message("  rows: ", n_rows)
if (n_rows == 0) stop("Nothing to plot.")
check_unmapped(rownames(dat$coef), metacyc_name_fn)

height <- min(30, max(6, 0.28 * n_rows + 3))
# Build the right-side category annotation once (row order = original matrix order;
# ComplexHeatmap will reorder it to match the main heatmap during draw()).
right_ann <- make_category_ann(rownames(dat$coef))

# Wider PDF to accommodate the longer English names + right annotation bar.
for (order_mode in c("padj", "clust")) {
  ht <- build_heatmap(
    dat,
    order_mode   = order_mode,
    title_prefix = sprintf(
      "PICRUSt2 pathway (padj<%.2f in any disease, %d features)",
      MAIN_FILTER_PADJ, n_rows
    ),
    row_label_fn   = metacyc_name_fn,
    thresholds     = thresholds,
    meta_indicator = "star"
  )
  # Attach right-side annotation.  The + operator creates a HeatmapList;
  # row ordering is inherited from ht (the first heatmap in the list).
  ht_combined <- ht + right_ann

  out <- file.path(HEATDIR,
                   sprintf("heatmap_picrust2_pathway_by_%s.pdf", order_mode))
  pdf(out, width = 15, height = height)   # +1 in to accommodate annotation bar
  ComplexHeatmap::draw(
    ht_combined,
    merge_legend           = TRUE,
    heatmap_legend_side    = "right",
    annotation_legend_side = "right"
  )
  dev.off()
  message("  wrote ", out)
}

message("\nDone. PICRUSt2 pathway heatmaps in ", HEATDIR)
message("To broaden filter: MAIN_FILTER_PADJ <- 0.25 (38 pathways) or <- 1 (all 377).")
