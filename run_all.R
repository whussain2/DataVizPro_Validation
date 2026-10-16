# ==============================================================================
# run_all.R - reproduce every number in the manuscript's Table 4 and Sections 3.5-3.7
# Usage (from the repository root):  Rscript run_all.R
# Total time: about 30 min on 2 cores (the EM ancestry engine takes most of it)
# ==============================================================================
scripts <- c("01_barley_spatial.R", "02_rice_qc.R", "03_ancestry.R", "04_kinship.R",
             "05_genotype_pca.R", "06_ld_heatmap.R", "07_ld_decay.R", "08_gwas.R",
             "09_phenotypes.R", "10_figures.R")
for (s in scripts) {
  message("\n==== ", s, " ====")
  t0 <- Sys.time()
  system2("Rscript", file.path("scripts", s))
  message(sprintf("%s finished in %.1f min", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo.txt"))
