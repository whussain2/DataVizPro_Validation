# ==============================================================================
# 03_ancestry.R - ancestry coefficients at K = 5 (Table 4, Section 3.6)
# DataVizPro engines: block-wise EM admixture (5 runs) and sNMF via LEA (K = 1..8)
# Reference: fastmixture 1.3 (ADMIXTURE model) - run python/references.py fastmixture
# ==============================================================================
source("scripts/00_setup.R")
prep10 <- readRDS(file.path(TMP_DIR, "prep10.rds"))

em_file <- file.path(TMP_DIR, "em_k5.rds")
if (!file.exists(em_file)) {          # about 20 min on 2 cores
  em <- ps_run_parametric_admixture(prep10, engine = "em", k_min = 5L, k_max = 5L,
                                    repetitions = 5L, iterations = 200L, seed = 42L)
  saveRDS(em, em_file)
}
em <- readRDS(em_file)
sn_file <- file.path(TMP_DIR, "snmf.rds")
if (!file.exists(sn_file)) {
  sn <- ps_run_snmf(prep10, k_min = 1L, k_max = 8L, repetitions = 5L, cpu = 1L, seed = 42L)
  saveRDS(sn, sn_file)
}
sn <- readRDS(sn_file)

fm_file <- file.path(TMP_DIR, "fm_k5.K5.s1.Q")
if (!file.exists(fm_file) && nzchar(Sys.which("fastmixture")) && nzchar(Sys.which("python3")))
  system2("python3", c("python/references.py", "fastmixture"))
stopifnot("fastmixture output missing: run python/references.py fastmixture" = file.exists(fm_file))
fm <- as.matrix(read.table(fm_file))

# Match clusters by optimal assignment (Hungarian algorithm; Kuhn 1955)
align <- function(ref, q) {
  s <- cor(ref, q); a <- clue::solve_LSAP(s - min(s), maximum = TRUE)
  q[, as.integer(a)]
}
q_em   <- em$q_by_k[["5"]]$best     # run with lowest masked cross-entropy (DataVizPro default)
q_snmf <- sn$q_by_k[["5"]]$best
a_em <- cmp(as.vector(q_em), as.vector(align(q_em, fm)))
a_sn <- cmp(as.vector(q_snmf), as.vector(align(q_snmf, fm)))
table4_row("Ancestry coefficients (EM engine)", "fastmixture 1.3 (ADMIXTURE model)", "Rice 44K, K = 5",
           a_em["n"], a_em["r"], a_em["max_abs_diff"], "10,000 markers; best of 5 runs")
table4_row("Ancestry coefficients (sNMF engine)", "fastmixture 1.3 (ADMIXTURE model)", "Rice 44K, K = 5",
           a_sn["n"], a_sn["r"], a_sn["max_abs_diff"], "10,000 markers; best of 5 runs")

# Ancestry groups (dominant coefficient >= 0.7) used by later scripts
rownames(q_em) <- prep10$ids
grp <- ifelse(apply(q_em, 1, max) >= 0.7, paste0("G", max.col(q_em)), "admixed")
saveRDS(grp, file.path(TMP_DIR, "ancestry_groups.rds"))
write.csv(data.frame(accession = prep10$ids, round(q_em, 4), group = grp),
          file.path(TAB_DIR, "case2_rice_ancestry_K5.csv"), row.names = FALSE)
write.csv(sn$summary, file.path(TAB_DIR, "case2_rice_snmf_cross_entropy.csv"), row.names = FALSE)
write.csv(em$runs, file.path(TAB_DIR, "case2_rice_em_runs.csv"), row.names = FALSE)
print(table(grp)); print(sn$summary); print(em$runs)
