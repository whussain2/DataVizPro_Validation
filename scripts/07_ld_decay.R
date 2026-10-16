# ==============================================================================
# 07_ld_decay.R - LD decay in the two largest ancestry groups (Table 4, Section 3.6)
# DataVizPro LD decay module versus an independent least-squares Hill-Weir fit
# ==============================================================================
source("scripts/00_setup.R")
G <- read_rice_geno(); map <- read_rice_map()
stopifnot(identical(colnames(G), map$snp))
grp <- readRDS(file.path(TMP_DIR, "ancestry_groups.rds"))
big <- names(sort(table(grp[grp != "admixed"]), decreasing = TRUE))[1:2]

rows <- list()
for (g in big) {
  ids <- names(grp)[grp == g]
  f <- ldd_filter_snps(G[ids, ], map$bp, map$chr, maf_min = 0.05, callrate_min = 0.8)
  ld <- do.call(rbind, lapply(1:12, function(ch) {
    k <- f$chr == ch
    ldd_compute_ld_pairs(f$X[, k], f$pos[k], max_dist = 1e6, max_snps = 6000L, pair_cap = 2e5, seed = ch)
  }))
  entry <- ldd_make_entry(ld, n_ind = length(ids))
  half_dvp <- ldd_half_decay(entry$fit, dist_max = 5e6)

  # Independent fit: minimise the residual sum of squares over log(rho), then solve for the half point
  n <- 2 * length(ids)
  hw <- function(d, rho) {
    C <- rho * d
    ((10 + C) / ((2 + C) * (11 + C))) * (1 + ((3 + C) * (12 + 12 * C + C^2)) / (n * (2 + C) * (11 + C)))
  }
  o <- optimize(function(lr) sum((ld$r2 - hw(ld$dist, exp(lr)))^2), c(log(1e-10), log(1e-1)), tol = 1e-10)
  rho <- exp(o$minimum); r0 <- hw(1, rho)
  half_ref <- uniroot(function(d) hw(d, rho) - r0 / 2, c(1, 5e7), tol = 1e-6)$root
  rows[[g]] <- data.frame(group = g, accessions = length(ids), snps = ncol(f$X), pairs = nrow(ld),
                          rho_datavizpro = entry$fit$rho, rho_reference = rho,
                          half_decay_kb_datavizpro = half_dvp / 1000, half_decay_kb_reference = half_ref / 1000)
  saveRDS(entry$binned, file.path(TMP_DIR, sprintf("ld_decay_binned_%s.rds", g)))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(TAB_DIR, "case2_rice_ld_decay.csv"), row.names = FALSE)
table4_row("LD decay fit (half-decay distance)", "Independent least-squares Hill-Weir fit", "Rice 44K, 2 ancestry groups",
           nrow(res), NA, max(abs(res$half_decay_kb_datavizpro - res$half_decay_kb_reference)), "difference in kb")
print(res)
