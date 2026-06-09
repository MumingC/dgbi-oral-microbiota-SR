# 00_setup.R — Install and load packages for DGBI oral microbiota analysis
# Run once to install, then source() at the top of each script.

cran_pkgs <- c("tidyverse", "vegan", "ape", "patchwork", "here",
               "broom", "ggpubr", "metafor", "remotes", "BiocManager")
bioc_pkgs <- c("phyloseq", "Maaslin2", "ANCOMBC", "microbiome")

installed <- rownames(installed.packages())
to_install_cran <- setdiff(cran_pkgs, installed)
if (length(to_install_cran)) install.packages(to_install_cran)

if (!"BiocManager" %in% installed) install.packages("BiocManager")
to_install_bioc <- setdiff(bioc_pkgs, installed)
if (length(to_install_bioc)) BiocManager::install(to_install_bioc, update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(vegan)
  library(ape)
  library(here)
})

# Project paths — always anchor to project root via here::here()
DATA_ROOT <- here::here("..")          # dgbi_exports/ (parent of R/)
RESULTS   <- here::here("results")
FIGURES   <- here::here("figures")
CACHE     <- here::here("cache")
dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE,   showWarnings = FALSE, recursive = TRUE)

message("Setup complete. DATA_ROOT = ", normalizePath(DATA_ROOT))
