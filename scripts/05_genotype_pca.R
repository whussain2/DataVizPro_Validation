# ==============================================================================
# 05_genotype_pca.R - genotypic principal components (Section 3.6)
# ==============================================================================
source("scripts/00_setup.R")
prep <- readRDS(file.path(TMP_DIR, "prep_all.rds"))
pc <- ps_run_pca(prep, components = 10L)
grp <- readRDS(file.path(TMP_DIR, "ancestry_groups.rds"))
scores <- pc$scores[, c("sample_id", paste0("PC", 1:5))]
scores$group <- grp[scores$sample_id]
write.csv(scores, file.path(TAB_DIR, "case2_rice_genotype_pcs.csv"), row.names = FALSE)
r2 <- function(pc) summary(lm(scores[[pc]] ~ scores$group))$r.squared
out <- data.frame(statistic = c("PC1 variance (%)", "PC2 variance (%)",
                                "R2 of PC1 on ancestry group", "R2 of PC2 on ancestry group"),
                  value = round(c(100 * pc$variance[1:2], r2("PC1"), r2("PC2")), 3))
write.csv(out, file.path(TAB_DIR, "case2_rice_genotype_pca.csv"), row.names = FALSE)
print(out)
