# 11_fig2_tally_combined.R
# Combined Figure 2 for Gut Microbes manuscript:
#   Panel (a): GERD-associated genera reported in >=2 included studies
#   Panel (b): IBS-associated genera reported in >=2 included studies
#
# Data source: MANUALLY-EXTRACTED matrices from each included paper's reported
# significant taxa (NOT the standardised MaAsLin2 re-analysis).
#   - GERD_genera_matrix.csv          (Source A for GERD)
#   - IBS_genera_difference_v2.csv    (Source A for IBS)
#
# Scoring rule (per the manuscript Methods):
#   "up"                 -> +1
#   "down"               -> -1
#   "up (mild)" etc.     -> +1   (severity qualifier ignored for tally)
#   "down (mild)" etc.   -> -1
#   "up/down"            -> +0.5 in the up bar AND -0.5 in the down bar
#                           (within-study discordant: same genus, different sites)
#
# Inclusion filter: genus appears in >= 2 different studies (>=2 non-empty cells).
#
# Output:
#   figures/main/Fig2_genera_tally_combined.pdf  (vector, journal-ready)
#   figures/main/Fig2_genera_tally_combined.tiff (300 dpi)
#   figures/main/Fig2_genera_tally_combined.png  (300 dpi, for Word review)
#
# Run from project root:
#   Rscript scripts/11_fig2_tally_combined.R

# Source project setup (loads tidyverse, ggplot2, patchwork, here, etc.)
source(here::here("00_setup.R"))

suppressPackageStartupMessages({
  library(stringr)
  library(patchwork)
})

# ---- Paths ------------------------------------------------------------------
# Project layout:
#   <parent>/dgbi_exports/R/scripts/11_fig2_tally_combined.R  (this script)
#   <parent>/Systemic review on oral microbiota of DGBI patients/  (data CSVs)
#
# here::here() resolves to dgbi_exports/R/ (anchored by 00_setup.R or .here file),
# so the data folder is two levels up + sibling of dgbi_exports/.
DATA_DIR <- here::here("..", "..",
                       "Systemic review on oral microbiota of DGBI patients")
OUT_DIR  <- file.path(FIGURES, "main")   # FIGURES is defined in 00_setup.R
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GERD_CSV <- file.path(DATA_DIR, "GERD_genera_matrix.csv")
IBS_CSV  <- file.path(DATA_DIR, "IBS_genera_difference_v2.csv")

# Sanity check — fail fast with a clear message if the data isn't where expected
stopifnot(
  "GERD_genera_matrix.csv not found — check DATA_DIR" = file.exists(GERD_CSV),
  "IBS_genera_difference_v2.csv not found — check DATA_DIR" = file.exists(IBS_CSV)
)
message("DATA_DIR resolved to: ", normalizePath(DATA_DIR))
message("OUT_DIR  resolved to: ", normalizePath(OUT_DIR))

# ---- Parse one cell into c(up_score, down_score) ----------------------------
parse_cell <- function(x) {
  if (is.na(x) || str_trim(x) == "") return(c(up = 0, down = 0))
  v <- tolower(str_trim(x))
  # within-study discordant: "up/down" or "down/up"
  if (str_detect(v, "up\\s*/\\s*down|down\\s*/\\s*up")) {
    return(c(up = 0.5, down = -0.5))
  }
  # standard up / down (severity qualifier in parentheses is ignored)
  if (str_detect(v, "^up"))   return(c(up = 1,  down = 0))
  if (str_detect(v, "^down")) return(c(up = 0,  down = -1))
  warning("Unrecognised cell value: '", x, "'")
  c(up = 0, down = 0)
}

# ---- Build long-format scores -----------------------------------------------
build_scores <- function(csv_path) {
  raw <- read_csv(csv_path, show_col_types = FALSE)
  # First column = Genus; remaining columns = study names
  study_cols <- setdiff(names(raw), names(raw)[1])
  long <- raw |>
    rename(genus = 1) |>
    pivot_longer(all_of(study_cols), names_to = "study", values_to = "raw") |>
    filter(!is.na(raw) & str_trim(raw) != "")

  # Apply scoring rule
  scored <- long |>
    rowwise() |>
    mutate(parsed = list(parse_cell(raw))) |>
    mutate(up = parsed[["up"]], down = parsed[["down"]]) |>
    ungroup() |>
    select(genus, study, up, down)
  scored
}

# ---- Aggregate per genus and filter to >=2 studies --------------------------
summarise_genera <- function(scored, min_studies = 2L) {
  per_genus <- scored |>
    group_by(genus) |>
    summarise(
      n_studies = n_distinct(study),
      up_sum    = sum(up),         # always >= 0 (positive contributions)
      down_sum  = sum(down),       # always <= 0 (down stored as negative)
      net_score = sum(up) + sum(down),  # signed net direction
      .groups   = "drop"
    ) |>
    mutate(magnitude = up_sum - down_sum) |>   # total bar length (down is negative, so this is up + |down|)
    filter(n_studies >= min_studies) |>
    # Primary: most positive net at top, most negative (e.g., Neisseria) at bottom.
    # Secondary (tie-breaker): larger total bar length first within the same net_score.
    # NOTE: down_sum is already negative, so net = up_sum + down_sum (NOT minus).
    arrange(desc(net_score), desc(magnitude))
  per_genus
}

# ---- Plot one panel ---------------------------------------------------------
tally_panel <- function(df, title, x_label = "Number of studies") {
  long <- df |>
    select(genus, up_sum, down_sum) |>
    pivot_longer(c(up_sum, down_sum), names_to = "direction", values_to = "value") |>
    mutate(
      direction = factor(direction, levels = c("up_sum", "down_sum"),
                         labels = c("Increased", "Decreased")),
      genus     = factor(genus, levels = rev(df$genus))   # preserve order
    )

  ggplot(long, aes(x = value, y = genus, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, color = "grey20") +
    scale_fill_manual(values = c(Increased = "#d4652a", Decreased = "#3970b0"),
                      name = NULL) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(title = title, x = x_label, y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      plot.title       = element_text(face = "bold", size = 11),
      axis.text.y      = element_text(face = "italic"),
      legend.position  = "bottom",
      legend.margin    = margin(0, 0, 0, 0),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
    )
}

# ---- Build the two panels ---------------------------------------------------
gerd_scores  <- build_scores(GERD_CSV)
ibs_scores   <- build_scores(IBS_CSV)

gerd_summary <- summarise_genera(gerd_scores, min_studies = 2L)
ibs_summary  <- summarise_genera(ibs_scores,  min_studies = 2L)

message(sprintf("GERD: %d genera in >=2 studies", nrow(gerd_summary)))
message(sprintf("IBS:  %d genera in >=2 studies", nrow(ibs_summary)))

p_a <- tally_panel(gerd_summary,
                   title = "(a) GERD-spectrum (n = 7 studies)")
p_b <- tally_panel(ibs_summary,
                   title = "(b) IBS (n = 4 studies)")

# Independent x-axes preserved (no shared scale).
# Heights weighted by number of genera so bars are comparable in thickness.
heights <- c(nrow(gerd_summary), nrow(ibs_summary))
fig2 <- p_a / p_b +
  plot_layout(heights = heights, guides = "collect") &
  theme(legend.position = "bottom")

# ---- Save -------------------------------------------------------------------
# Total figure height scaled to total number of genera (0.25 in per genus + padding)
total_h <- 0.25 * sum(heights) + 2
ggsave(file.path(OUT_DIR, "Fig2_genera_tally_combined.pdf"),
       fig2, width = 7, height = total_h, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Fig2_genera_tally_combined.tiff"),
       fig2, width = 7, height = total_h, dpi = 300, compression = "lzw")
ggsave(file.path(OUT_DIR, "Fig2_genera_tally_combined.png"),
       fig2, width = 7, height = total_h, dpi = 300)

message("\nWrote 3 files to: ", normalizePath(OUT_DIR))
message("  Fig2_genera_tally_combined.pdf  (vector, for submission)")
message("  Fig2_genera_tally_combined.tiff (300 dpi, for submission)")
message("  Fig2_genera_tally_combined.png  (300 dpi, for Word review)")
