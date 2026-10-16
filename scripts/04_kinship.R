# ==============================================================================
# 04_kinship.R - VanRaden (2008) genomic relationship matrix (Table 4, Section 3.6)
# DataVizPro: gh_kinship(); references: the published formula and rrBLUP::A.mat
# ==============================================================================
source("scripts/00_setup.R")
prep <- readRDS(file.path(TMP_DIR, "prep_all.rds"))
M <- prep$geno
K_dvp <- gh_kinship(M, "vanraden")

# Reference 1: VanRaden method 1 written out, with mean imputation
X <- M
for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
p <- colMeans(X) / 2
Z <- sweep(X, 2, 2 * p)
G_ref <- tcrossprod(Z) / (2 * sum(p * (1 - p)))
# Reference 2: rrBLUP (genotypes coded -1/0/1)
A <- rrBLUP::A.mat(M - 1, min.MAF = 0, max.missing = 1, impute.method = "mean", shrink = FALSE)

off <- upper.tri(K_dvp)
a1 <- cmp(K_dvp[off], G_ref[off]); a2 <- cmp(K_dvp[off], A[off])
table4_row("Genomic relationship matrix", "VanRaden formula (base R); rrBLUP 4.6", "Rice 44K",
           a1["n"], min(a1["r"], a2["r"]), max(a1["max_abs_diff"], a2["max_abs_diff"]),
           sprintf("%d markers; off-diagonal elements", ncol(M)))

grp <- readRDS(file.path(TMP_DIR, "ancestry_groups.rds"))[rownames(K_dvp)]
pure <- grp != "admixed"
same <- outer(grp, grp, "==") & outer(pure, pure, "&")
diff <- !outer(grp, grp, "==") & outer(pure, pure, "&")
write.csv(data.frame(statistic = c("Mean relationship within ancestry groups", "Mean relationship between ancestry groups"),
                     value = round(c(mean(K_dvp[off & same]), mean(K_dvp[off & diff])), 3)),
          file.path(TAB_DIR, "case2_rice_kinship.csv"), row.names = FALSE)
saveRDS(K_dvp, file.path(TMP_DIR, "grm.rds"))
print(rbind(a1, a2))
