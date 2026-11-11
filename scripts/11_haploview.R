# ==============================================================================
# 11_haploview.R - HaploView module versus Haploview 4.2 (Supplemental Table 4)
#
# Step 1 (this script, no Haploview needed): writes a Haploview input file pair
#   for one chromosome of the public rice 44K panel and computes D', r2, LOD,
#   the D' interval and the three block definitions with the DataVizPro engine
#   (dvp_engine/haploview.R).
#
# Step 2 (Haploview 4.2, Java, run once on your computer):
#   java -jar Haploview.jar -nogui -pedfile results/haploview/rice_chr4.ped \
#        -info results/haploview/rice_chr4.info -dprime -blockoutput ALL \
#        -maxdistance 500 -minMAF 0 -hwcutoff 0 -missingCutoff 1 \
#        -out results/haploview/haploview42
#   This writes haploview42.LD, haploview42.GABRIELblocks,
#   haploview42.SPINEblocks and haploview42.4GAMblocks.
#
# Step 3 (this script again): if the Haploview files are present, the two sets
#   of results are compared pair by pair and block by block, and the rows are
#   added to results/tables/table4_agreement.csv (the rows of the paper's
#   747-line comparison are kept).
#
# Usage:  Rscript scripts/11_haploview.R [chromosome]   (default chromosome 4)
# ==============================================================================

source("scripts/00_setup.R")
args <- commandArgs(trailingOnly = TRUE)
CHR  <- if (length(args)) as.integer(args[1]) else 4L
HV_DIR <- file.path("results", "haploview")
dir.create(HV_DIR, showWarnings = FALSE, recursive = TRUE)
base <- file.path(HV_DIR, sprintf("rice_chr%d", CHR))

# ---- genotypes of one chromosome (0/1/2), same QC as the other analyses ------
G   <- read_rice_geno()
map <- read.table(file.path(DATA_DIR, "rice", "sativas413.map"),
                  col.names = c("chr", "snp", "cM", "bp"), stringsAsFactors = FALSE)
map <- map[map$snp %in% colnames(G) & map$chr == CHR, ]
map <- map[order(map$bp), ]
G   <- G[, map$snp, drop = FALSE]
call_rate <- colMeans(!is.na(G))
p   <- colMeans(G, na.rm = TRUE) / 2
maf <- pmin(p, 1 - p)
ok  <- call_rate >= 0.8 & maf >= 0.05
G   <- G[, ok, drop = FALSE]; map <- map[ok, ]
cat(sprintf("Chromosome %d: %d lines x %d markers after QC\n", CHR, nrow(G), ncol(G)))

# ---- Haploview input: linkage-format pedigree file and marker info file ------
# 0 -> "1 1", 1 -> "1 2", 2 -> "2 2", missing -> "0 0"
code <- c("1 1", "1 2", "2 2")
geno_txt <- apply(G, 1, function(x) {
  y <- ifelse(is.na(x), "0 0", code[x + 1])
  paste(y, collapse = " ")
})
ids <- rownames(G)
ped <- paste(ids, ids, 0, 0, 0, 0, geno_txt)
writeLines(ped, paste0(base, ".ped"))
write.table(map[, c("snp", "bp")], paste0(base, ".info"),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# ---- DataVizPro HaploView engine ----------------------------------------------
eng <- hv_ld_engine(G, pos_bp = map$bp, max_dist_kb = 500)
st  <- eng$st
dvp_blocks <- list(
  GABRIEL = hv_blocks_gabriel(st, eng$pos, eng$maf),
  SPINE   = hv_blocks_spine(st, eng$pos, max_dist_kb = 500),
  "4GAM"  = hv_blocks_4gamete(st, eng$pos)
)
out <- data.frame(L1 = eng$markers[st$i], L2 = eng$markers[st$j],
                  Dp = st$Dp, LOD = st$LOD, r2 = st$r2, CIlow = st$CIlow, CIhi = st$CIhi)
write.csv(out, paste0(base, "_datavizpro_LD.csv"), row.names = FALSE)
for (b in names(dvp_blocks)) {
  bl <- vapply(dvp_blocks[[b]], function(m) paste(eng$markers[m], collapse = " "), "")
  writeLines(bl, paste0(base, "_datavizpro_", b, "blocks.txt"))
}
cat(sprintf("DataVizPro: %d pairs within 500 kb; blocks: Gabriel %d, spine %d, four-gamete %d\n",
            nrow(st), length(dvp_blocks$GABRIEL), length(dvp_blocks$SPINE), length(dvp_blocks[["4GAM"]])))

# ---- comparison with Haploview 4.2, when its output is present -----------------
hv_ld_file <- file.path(HV_DIR, "haploview42.LD")
if (!file.exists(hv_ld_file)) {
  cat("Haploview 4.2 output not found in", HV_DIR, "- run step 2 (see the header), then rerun this script.\n")
  quit(save = "no")
}

hv <- read.table(hv_ld_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
names(hv)[1:7] <- c("L1", "L2", "Dp", "LOD", "r2", "CIlow", "CIhi")
key_dvp <- paste(out$L1, out$L2)
key_hv  <- paste(hv$L1, hv$L2)
m <- match(key_hv, key_dvp)
both <- !is.na(m)
cat(sprintf("Pairs in Haploview output: %d; matched: %d\n", nrow(hv), sum(both)))

# Haploview prints 2 decimals (D', r2, interval) and 2 decimals for LOD:
# agreement "within printed rounding" means |difference| <= 0.005 (+ float slack)
tol <- 0.005 + 1e-9
for (v in c("Dp", "r2", "LOD", "CIlow", "CIhi")) {
  a <- out[[v]][m[both]]; b <- hv[[v]][both]
  d <- abs(a - b)
  cat(sprintf("%-5s max |DataVizPro - Haploview printed value| = %.4f; pairs outside rounding = %d\n",
              v, max(d, na.rm = TRUE), sum(abs(a - b) > tol, na.rm = TRUE)))
}
cmp_dp <- cmp(out$Dp[m[both]], hv$Dp[both])
table4_row(sprintf("HaploView D', r2, LOD and D' interval (rice 44K chr %d)", CHR), "Haploview 4.2",
           sprintf("Rice 44K, chr %d (%d lines x %d markers)", CHR, nrow(G), ncol(G)),
           sum(both), cmp_dp["r"], cmp_dp["max_abs_diff"],
           "D' shown; all statistics compared after rounding to Haploview's printed precision")

# Haploview block files: lines starting with "BLOCK" list marker numbers (1-based, .info order)
read_hv_blocks <- function(f) {
  if (!file.exists(f)) return(NULL)
  x <- grep("^BLOCK", readLines(f), value = TRUE)
  lapply(strsplit(sub("^BLOCK [0-9]+\\.\\s+MARKERS:\\s*", "", x), "\\s+"),
         function(z) map$snp[as.integer(sub("!$", "", z))])
}
for (b in names(dvp_blocks)) {
  hb <- read_hv_blocks(file.path(HV_DIR, paste0("haploview42.", b, "blocks")))
  if (is.null(hb)) next
  db <- lapply(dvp_blocks[[b]], function(mm) eng$markers[mm])
  same <- sum(vapply(hb, function(h) any(vapply(db, identical, TRUE, h)), TRUE))
  cat(sprintf("%-8s blocks: Haploview %d, DataVizPro %d, identical %d\n", b, length(hb), length(db), same))
  table4_row(sprintf("HaploView haplotype blocks, %s (rice 44K chr %d)", b, CHR), "Haploview 4.2 (default settings)",
             sprintf("Rice 44K, chr %d", CHR), length(hb), NA, NA,
             sprintf("identical: %d of %d (DataVizPro found %d)", same, length(hb), length(db)))
}
