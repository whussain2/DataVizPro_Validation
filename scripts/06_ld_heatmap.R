# ==============================================================================
# 06_ld_heatmap.R - pairwise LD r2 and D' (Supplemental Table 4)
# DataVizPro LD heatmap module (ld_compute) versus snpStats (EM haplotype estimates)
# and the squared genotype correlation (the PLINK --r2 definition)
# 500 consecutive QC'd markers on chromosome 1
# ==============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages(library(snpStats))
prep <- readRDS(file.path(TMP_DIR, "prep_all.rds"))
map <- read_rice_map()
chr <- map$chr[match(colnames(prep$geno), map$snp)]
M <- prep$geno[, which(chr == 1)[1:500]]

r2_dvp <- ld_compute(M, "r2")
dp_dvp <- ld_compute(M, "dprime")
sm <- new("SnpMatrix", matrix(as.raw(ifelse(is.na(M), 0, M + 1)), nrow(M), dimnames = dimnames(M)))
ref <- ld(sm, stats = c("R.squared", "D.prime"), depth = ncol(M) - 1)
r2_em <- as.matrix(ref$R.squared); dp_em <- as.matrix(ref$D.prime)
r2_corr <- suppressWarnings(cor(M, use = "pairwise.complete.obs")^2)
ut <- upper.tri(r2_em)
a1 <- cmp(r2_dvp[ut], r2_em[ut]); a2 <- cmp(r2_dvp[ut], r2_corr[ut]); a3 <- cmp(dp_dvp[ut], dp_em[ut])
table4_row("LD r2", "snpStats 1.62 (EM r2)", "Rice 44K, chr 1 (500 SNPs)", a1["n"], a1["r"], a1["max_abs_diff"],
           "differences only for pairs with heterozygous or missing calls")
table4_row("LD r2", "Squared genotype correlation (PLINK --r2)", "Rice 44K, chr 1 (500 SNPs)", a2["n"], a2["r"], a2["max_abs_diff"], "")
table4_row("LD D'", "snpStats 1.62 (EM D')", "Rice 44K, chr 1 (500 SNPs)", a3["n"], a3["r"], a3["max_abs_diff"], "")
print(rbind(a1, a2, a3))
