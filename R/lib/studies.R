# lib/studies.R — Study registry: disease, path, region, sample counts.
# One row per study. Used by all downstream scripts to iterate consistently.

studies <- tibble::tribble(
  ~study,             ~disease, ~dir,                       ~region, ~n_samples, ~layout,
  "Ziganshina_2020",  "GERD",   "GERD/Ziganshina_2020",     "V3-V4", 26L,        "PE",
  "Qian_2023",        "GERD",   "GERD/Qian_2023",           "V3-V4", 60L,        "PE",
  "Li_2025",          "IBS-D",  "IBS/Li_2025",              "V3-V4", 35L,        "PE",
  "Zheng_2024",       "LPRD",   "GERD/Zheng_2024",          "V3-V4", 182L,       "PE",
  "Hao_2022",         "GERD",   "GERD/Hao_2022",            "V3-V4", 24L,        "SE",
  "Kawar_2021",       "GERD",   "GERD/Kawar_2021",          "V1-V3", 128L,       "SE",
  "Tang_2023",        "IBS-D",  "IBS/Tang_2023",            "V4-V5", 62L,        "PE"
)
