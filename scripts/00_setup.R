# ==============================================================================
# 00_setup.R - paths, packages and the DataVizPro engine
# Run every script from the repository root, e.g.  Rscript scripts/01_barley_spatial.R
# ==============================================================================

ENGINE_DIR <- "dvp_engine"
DATA_DIR   <- "data"
TAB_DIR    <- file.path("results", "tables")
FIG_DIR    <- file.path("results", "figures")
TMP_DIR    <- file.path("results", "intermediate")   # large objects, not committed
for (d in c(TAB_DIR, FIG_DIR, TMP_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# The unchanged DataVizPro v1.1.1 computation functions
source(file.path(ENGINE_DIR, "load_engine.R"))

# Agreement summary used throughout Table 4: Pearson r and maximum absolute difference
cmp <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  c(n = sum(ok), r = cor(a[ok], b[ok]), max_abs_diff = max(abs(a[ok] - b[ok])))
}

# Append one row to the Table 4 results file
table4_row <- function(analysis, reference, dataset, n, r, max_abs_diff, note = "") {
  f <- file.path(TAB_DIR, "table4_agreement.csv")
  row <- data.frame(analysis, reference, dataset, n, r = signif(r, 6),
                    max_abs_diff = signif(max_abs_diff, 4), note)
  old <- if (file.exists(f)) read.csv(f) else NULL
  if (!is.null(old)) old <- old[!(old$analysis == analysis & old$reference == reference), ]
  write.csv(rbind(old, row), f, row.names = FALSE)
}

# Rice genotypes (413 accessions x 36,901 SNPs, 0/1/2) and the MSU6 marker map
read_rice_geno <- function() {
  f <- file.path(TMP_DIR, "G_raw.rds")
  if (file.exists(f)) return(readRDS(f))
  path <- file.path(DATA_DIR, "rice", "Geno44k.txt")
  hdr <- scan(path, what = "", nlines = 1, sep = "\t", quiet = TRUE)
  g <- data.table::fread(path, header = FALSE, skip = 1, na.strings = "NA")
  G <- as.matrix(g[, -1]); storage.mode(G) <- "double"
  dimnames(G) <- list(g[[1]], hdr)
  saveRDS(G, f); G
}
read_rice_map <- function() {
  read.table(file.path(DATA_DIR, "rice", "sativas413.map"), col.names = c("chr", "snp", "cm", "bp"))
}
