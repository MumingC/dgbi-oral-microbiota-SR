# 01_build_phyloseq.R — Build SILVA + eHOMD phyloseq objects per study, cache as .rds.

source(here::here("00_setup.R"))
source(here::here("lib", "studies.R"))
source(here::here("lib", "build_phyloseq.R"))
source(here::here("lib", "filter_phyloseq.R"))

ps_list <- list()
for (i in seq_len(nrow(studies))) {
  s   <- studies$study[i]
  dir <- file.path(DATA_ROOT, studies$dir[i])
  for (db in c("silva", "ehomd")) {
    key <- paste(s, db, sep = "__")
    message("Building: ", key)
    ps_list[[key]] <- tryCatch({
      ps <- build_phyloseq(dir, db = db)
      n0 <- phyloseq::nsamples(ps)
      ps <- filter_phyloseq(ps, min_reads = 2000)
      message("  filtered: ", n0, " → ", phyloseq::nsamples(ps), " samples, ",
              phyloseq::ntaxa(ps), " taxa")
      ps
    }, error = function(e) { warning(key, ": ", conditionMessage(e)); NULL })
  }
}

saveRDS(ps_list, file.path(CACHE, "ps_list.rds"))
message("Cached ", length(ps_list), " phyloseq objects to ", file.path(CACHE, "ps_list.rds"))
