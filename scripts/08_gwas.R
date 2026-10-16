# ==============================================================================
# 08_gwas.R - GWAS summary statistics visualised in Fig. 6F (Section 3.6)
# Computed externally, as in the manuscript: rrBLUP mixed model (P3D), genomic
# relationship matrix + 3 principal components, MAF >= 0.05, Bonferroni alpha = 0.05
# ==============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages(library(rrBLUP))
prep <- readRDS(file.path(TMP_DIR, "prep_all.rds"))
map <- read_rice_map()
M <- prep$geno - 1
mm <- map[match(colnames(M), map$snp), ]
geno <- data.frame(marker = mm$snp, chrom = mm$chr, pos = mm$bp, t(M), check.names = FALSE)
ph <- read.csv(file.path(DATA_DIR, "rice", "rice_phenotypes.csv"), row.names = 1)

loci_all <- list(); summ <- list()
for (tr in c("Seed.length", "Flowering.time.at.Arkansas")) {
  pheno <- data.frame(gid = rownames(ph), y = ph[[tr]]); names(pheno)[2] <- tr
  pheno <- pheno[pheno$gid %in% rownames(M) & is.finite(pheno[[tr]]), ]
  g <- GWAS(pheno, geno, n.PC = 3, min.MAF = 0.05, P3D = TRUE, plot = FALSE)
  names(g)[4] <- "logP"
  thr <- -log10(0.05 / nrow(g))
  sig <- g[g$logP > thr, ]; sig <- sig[order(-sig$logP), ]
  loci <- sig[0, ]   # merge significant SNPs within 1 Mb into loci
  for (i in seq_len(nrow(sig)))
    if (!any(loci$chrom == sig$chrom[i] & abs(loci$pos - sig$pos[i]) < 1e6)) loci <- rbind(loci, sig[i, ])
  lambda <- median(qchisq(10^-g$logP, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
  write.csv(data.frame(SNP = g$marker, CHR = g$chrom, BP = g$pos, P = 10^-g$logP),
            file.path(TMP_DIR, sprintf("gwas_%s.csv", tr)), row.names = FALSE)
  loci_all[[tr]] <- data.frame(trait = tr, loci, row.names = NULL)
  summ[[tr]] <- data.frame(trait = tr, n = nrow(pheno), threshold_logP = round(thr, 2),
                           significant_snps = nrow(sig), chromosomes = length(unique(sig$chrom)),
                           loci_1Mb = nrow(loci), lambda = round(lambda, 3))
}
write.csv(do.call(rbind, summ), file.path(TAB_DIR, "case2_rice_gwas_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, loci_all), file.path(TAB_DIR, "case2_rice_gwas_loci.csv"), row.names = FALSE)
print(do.call(rbind, summ))
