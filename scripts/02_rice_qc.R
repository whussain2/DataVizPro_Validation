# ==============================================================================
# 02_rice_qc.R - genotype quality control with the DataVizPro population-structure
# preparation step (ps_prepare_matrix): call rate >= 0.8, MAF >= 0.05
# ==============================================================================
source("scripts/00_setup.R")
G <- read_rice_geno()
cat("Raw matrix:", nrow(G), "accessions x", ncol(G), "SNPs; missing", round(100 * mean(is.na(G)), 2), "%\n")

prep_all <- ps_prepare_matrix(G, ids = rownames(G), min_sample_call = 0.8,
                              min_marker_call = 0.8, min_maf = 0.05)
# The EM / variational-Bayes engines are limited to 20,000 markers: DataVizPro's own
# random marker cap (seed 42) selects 10,000 for the ancestry comparisons
prep10 <- ps_prepare_matrix(G, ids = rownames(G), min_sample_call = 0.8, min_marker_call = 0.8,
                            min_maf = 0.05, max_markers = 10000, seed = 42L)
saveRDS(prep_all, file.path(TMP_DIR, "prep_all.rds"))
saveRDS(prep10,   file.path(TMP_DIR, "prep10.rds"))
# Plain-text copy of the 10,000-marker matrix for the Python references (fastmixture)
write.table(prep10$geno, file.path(TMP_DIR, "prep10.tsv"), sep = "\t", quote = FALSE, na = "-9")

qc <- data.frame(step = c("Raw", "After QC", "Ancestry subset"),
                 accessions = c(nrow(G), nrow(prep_all$geno), nrow(prep10$geno)),
                 snps = c(ncol(G), ncol(prep_all$geno), ncol(prep10$geno)))
write.csv(qc, file.path(TAB_DIR, "case2_rice_qc.csv"), row.names = FALSE)
print(qc)
