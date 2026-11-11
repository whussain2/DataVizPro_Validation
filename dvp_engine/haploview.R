# ==============================================================================
# FILE: modules/genomicvisuals/haploview/hv_stats.R
# PURPOSE (AWS-012): HaploView statistics computed the way Haploview 4.2 does.
#   * hv_pair_stats()    two-locus EM haplotype frequencies per marker pair ->
#                        D', r2, LOD and the 90% likelihood interval of D'
#                        (Haploview's "95% confidence bounds": 5% in each tail)
#   * hv_ld_engine()     all pairs of a genotype matrix, with a maximum
#                        distance for the interval (Haploview -maxdistance)
#   * hv_blocks_gabriel() / hv_blocks_spine() / hv_blocks_4gamete()
#                        Haploview's three block definitions with its defaults
#   * hv_color_standard() Haploview's standard D'/LOD colour scheme
# Validated against Haploview 4.2 on 747 lines x 1,279 markers (chr 4):
# D', r2, LOD and the interval agree within Haploview's printed rounding for
# 817,281 pairs; Gabriel 184/184, solid spine 361/361 and four-gamete 329/329
# blocks are identical (max distance 500 kb).
# Genotypes are 0/1/2 dosages; anything else is treated as missing.
# ASCII only.
# ==============================================================================

# ---- pair statistics ---------------------------------------------------------
# cnt: matrix (pairs x 9) of two-locus genotype counts, columns
#      g1g2 = 00 01 02 10 11 12 20 21 22 (dosage of marker 1, marker 2).
# ci:  logical vector, TRUE where the D' interval is wanted.
.hv_stats_from_counts <- function(cnt, ci, max_iter = 1000L, tol = 1e-10) {
  storage.mode(cnt) <- "double"
  np <- nrow(cnt)
  k00 <- 2 * cnt[, 1] + cnt[, 2] + cnt[, 4]
  k01 <- 2 * cnt[, 3] + cnt[, 2] + cnt[, 6]
  k10 <- 2 * cnt[, 7] + cnt[, 4] + cnt[, 8]
  k11 <- 2 * cnt[, 9] + cnt[, 8] + cnt[, 6]
  dh  <- cnt[, 5]                                   # double heterozygotes (phase unknown)
  tot <- k00 + k01 + k10 + k11 + 2 * dh
  tot[tot <= 0] <- NA_real_
  f00 <- (k00 + .5 * dh) / tot; f01 <- (k01 + .5 * dh) / tot
  f10 <- (k10 + .5 * dh) / tot; f11 <- (k11 + .5 * dh) / tot
  # EM only where phase is ambiguous; pairs leave the active set once converged
  act <- which(dh > 0 & is.finite(tot))
  it <- 0L
  while (length(act) && it < max_iter) {
    it <- it + 1L
    a00 <- f00[act]; a01 <- f01[act]; a10 <- f10[act]; a11 <- f11[act]
    den <- a00 * a11 + a01 * a10
    cp  <- ifelse(den > 0, a00 * a11 / den, .5)
    d_  <- dh[act]; t_ <- tot[act]
    g00 <- (k00[act] + cp * d_) / t_;       g11 <- (k11[act] + cp * d_) / t_
    g01 <- (k01[act] + (1 - cp) * d_) / t_; g10 <- (k10[act] + (1 - cp) * d_) / t_
    ch  <- pmax(abs(g00 - a00), abs(g11 - a11), abs(g01 - a01), abs(g10 - a10))
    f00[act] <- g00; f01[act] <- g01; f10[act] <- g10; f11[act] <- g11
    act <- act[is.finite(ch) & ch >= tol]
  }
  # orient so that D >= 0 (Haploview swaps the alleles of marker 2)
  D   <- f00 * f11 - f01 * f10
  neg <- which(D < 0)
  if (length(neg)) {
    t <- f00[neg]; f00[neg] <- f01[neg]; f01[neg] <- t
    t <- f10[neg]; f10[neg] <- f11[neg]; f11[neg] <- t
    t <- k00[neg]; k00[neg] <- k01[neg]; k01[neg] <- t
    t <- k10[neg]; k10[neg] <- k11[neg]; k11[neg] <- t
  }
  D  <- abs(D)
  p1 <- f00 + f01; q1 <- 1 - p1; p2 <- f00 + f10; q2 <- 1 - p2
  Dmax <- pmin(p1 * q2, q1 * p2)
  vv   <- p1 * q1 * p2 * q2
  dp <- ifelse(is.finite(Dmax) & Dmax > 0, D / Dmax, 0)
  r2 <- ifelse(is.finite(vv) & vv > 0, D^2 / vv, 0)
  sl <- function(n, x) ifelse(n == 0, 0, n * log(pmax(x, 1e-300)))
  LL <- function(a, b, c, e, ix = TRUE)
    (sl(k00[ix], a) + sl(k01[ix], b) + sl(k10[ix], c) + sl(k11[ix], e) + sl(dh[ix], a * e + b * c)) / log(10)
  ll1 <- LL(f00, f01, f10, f11)
  ll0 <- LL(p1 * p2, p1 * q2, q1 * p2, q1 * q2)
  lod <- ll1 - ll0
  lo <- rep(NA_real_, np); hi <- rep(NA_real_, np)
  ix <- which(ci & is.finite(tot))
  if (length(ix)) {
    # likelihood of D' = 0, 0.01, ..., 1 (Haploview's 101-point surface)
    S <- matrix(0, length(ix), 101L)
    Dm <- Dmax[ix]; P1 <- p1[ix]; P2 <- p2[ix]; Q1 <- q1[ix]
    for (s in 0:100) {
      a <- (s / 100) * Dm + P1 * P2; b <- P1 - a; c <- P2 - a; e <- Q1 - c
      if (s == 100) { a <- pmax(a, 1e-10); b <- pmax(b, 1e-10); c <- pmax(c, 1e-10); e <- pmax(e, 1e-10) }
      S[, s + 1L] <- LL(a, b, c, e, ix)
    }
    Pm <- 10^(S - ll1[ix]); Pm[!is.finite(Pm)] <- 0
    tp <- rowSums(Pm)
    lowi <- rep(NA_integer_, length(ix)); cs <- numeric(length(ix))
    for (s in 1:101) { cs <- cs + Pm[, s]; hit <- is.na(lowi) & cs > 0.05 * tp; lowi[hit] <- s }
    highi <- rep(NA_integer_, length(ix)); cs <- numeric(length(ix))
    for (s in 101:1) { cs <- cs + Pm[, s]; hit <- is.na(highi) & cs > 0.05 * tp; highi[hit] <- s }
    lowi[is.na(lowi)] <- 1L; highi[is.na(highi)] <- 101L
    lo[ix] <- (lowi - 2L) / 100                     # first index past 5%, minus one step
    # (as in Haploview this can print -0.01 when D' = 0; it never affects blocks)
    hi[ix] <- pmin(highi, 100L) / 100               # (0-based index + 1), capped at 1
  }
  out <- data.frame(Dp = dp, r2 = r2, LOD = lod, CIlow = lo, CIhi = hi,
                    f00 = f00, f01 = f01, f10 = f10, f11 = f11)
  # a pair without information (one marker monomorphic among the lines typed at
  # both): Haploview reports D' 1, LOD 0, interval 0-0; block rules skip it
  degen <- !(is.finite(vv) & vv > 0)
  out$Dp[degen] <- 1; out$r2[degen] <- 0; out$LOD[degen] <- 0
  out$CIlow[degen & ci] <- 0; out$CIhi[degen & ci] <- 0
  out
}

# Two-locus genotype counts for pairs (i, j): column g1g2 of the result is the
# number of individuals with dosage g1 at marker i and g2 at marker j.
# ij / ji are the linear indices of (i, j) and (j, i) in a k x k matrix.
# Six cross-products suffice: count(b, a) at (i, j) = count(a, b) at (j, i).
.hv_pair_counts <- function(G, ij, ji) {
  I <- lapply(0:2, function(v) { x <- (G == v); x[is.na(x)] <- FALSE; storage.mode(x) <- "double"; x })
  cnt <- matrix(0L, length(ij), 9L, dimnames = list(NULL, c("00","01","02","10","11","12","20","21","22")))
  for (a in 0:2) for (b in a:2) {
    C <- if (a == b) crossprod(I[[a + 1L]]) else crossprod(I[[a + 1L]], I[[b + 1L]])
    cnt[, a * 3L + b + 1L] <- as.integer(C[ij])
    if (a != b) cnt[, b * 3L + a + 1L] <- as.integer(C[ji])
    rm(C)
  }
  cnt
}

#' Haploview statistics for marker pairs of a 0/1/2 genotype matrix.
#' @param M      individuals x markers dosage matrix
#' @param pairs  two-column integer matrix (i < j); default all pairs
#' @param ci     logical per pair (or a single value): compute the D' interval
#' @return data.frame(i, j, Dp, r2, LOD, CIlow, CIhi, f00, f01, f10, f11)
hv_pair_stats <- function(M, pairs = NULL, ci = TRUE) {
  G <- round(as.matrix(M)); storage.mode(G) <- "double"
  G[!is.finite(G) | G < 0 | G > 2] <- NA
  m <- ncol(G)
  if (is.null(pairs)) pairs <- which(upper.tri(matrix(FALSE, m, m)), arr.ind = TRUE)
  i <- pairs[, 1]; j <- pairs[, 2]
  ci <- rep_len(as.logical(ci), length(i))
  cnt <- .hv_pair_counts(G, (j - 1L) * m + i, (i - 1L) * m + j)
  st <- .hv_stats_from_counts(cnt, ci)
  cbind(data.frame(i = i, j = j), st)
}

# ---- all pairs of a region ---------------------------------------------------
#' @param geno         individuals x markers dosage matrix (monomorphic markers
#'                     are dropped first, as LD is undefined for them)
#' @param pos_bp       marker positions (bp) or NULL; incomplete -> marker order
#' @param chr          chromosome per marker or NULL (pairs on different
#'                     chromosomes are treated as unlinked: no interval)
#' @param max_dist_kb  pairs further apart get no D' interval and are not used
#'                     for blocks (Haploview -maxdistance, default 500; 0 = none)
#' @param progress     optional function(fraction) for a progress bar
#' @return list(keep, dp, r2, lod (k x k matrices), st (pairs within the
#'              distance, with interval and haplotype frequencies), maf, pos,
#'              positions_used ("bp" | "order"), max_dist_kb, n_pairs)
hv_ld_engine <- function(geno, pos_bp = NULL, chr = NULL, max_dist_kb = 500, progress = NULL) {
  G <- as.matrix(geno); storage.mode(G) <- "double"
  G <- round(G); G[!is.finite(G) | G < 0 | G > 2] <- NA
  keep <- apply(G, 2, function(x) { xx <- x[!is.na(x)]; length(xx) >= 2 && length(unique(xx)) > 1 })
  if (sum(keep) < 2)
    stop("At least two polymorphic markers are required to compute LD. ",
         "After dropping monomorphic markers, fewer than two remain.")
  G <- G[, keep, drop = FALSE]; k <- ncol(G); mk <- colnames(G)
  pos <- if (!is.null(pos_bp) && length(pos_bp) == length(keep)) suppressWarnings(as.numeric(pos_bp))[keep] else NULL
  used <- "bp"
  if (is.null(pos) || !all(is.finite(pos))) { pos <- as.numeric(seq_len(k)); used <- "order" }
  ch <- if (!is.null(chr) && length(chr) == length(keep)) as.character(chr)[keep] else rep("1", k)
  ch[is.na(ch) | !nzchar(ch)] <- "1"
  md <- suppressWarnings(as.numeric(max_dist_kb))
  if (length(md) != 1 || !is.finite(md) || md <= 0) md <- Inf
  p <- colMeans(G, na.rm = TRUE) / 2; maf <- pmin(p, 1 - p)

  ut <- which(upper.tri(matrix(FALSE, k, k)))           # column-major, i < j
  ii <- ((ut - 1L) %% k) + 1L; jj <- ((ut - 1L) %/% k) + 1L
  cnt <- .hv_pair_counts(G, ut, (ii - 1L) * k + jj)
  within <- ch[ii] == ch[jj] & abs(pos[jj] - pos[ii]) <= md * 1000
  dp <- matrix(0, k, k, dimnames = list(mk, mk)); r2 <- dp; lod <- dp
  st_parts <- list()
  chunk <- 50000L
  idx <- split(seq_along(ut), ceiling(seq_along(ut) / chunk))
  for (q in seq_along(idx)) {
    ix <- idx[[q]]
    s <- .hv_stats_from_counts(cnt[ix, , drop = FALSE], within[ix])
    dp[ut[ix]] <- s$Dp; r2[ut[ix]] <- s$r2; lod[ut[ix]] <- s$LOD
    w <- which(within[ix])
    if (length(w)) st_parts[[length(st_parts) + 1L]] <-
      cbind(data.frame(i = ii[ix][w], j = jj[ix][w]), s[w, , drop = FALSE])
    if (is.function(progress)) progress(q / length(idx))
  }
  rm(cnt)
  sym <- function(X) { X[lower.tri(X)] <- t(X)[lower.tri(X)]; diag(X) <- 1; X }
  dp <- sym(dp); r2 <- sym(r2); lod <- sym(lod); diag(lod) <- NA
  st <- if (length(st_parts)) do.call(rbind, st_parts) else
    data.frame(i = integer(0), j = integer(0), Dp = numeric(0), r2 = numeric(0), LOD = numeric(0),
               CIlow = numeric(0), CIhi = numeric(0), f00 = numeric(0), f01 = numeric(0), f10 = numeric(0), f11 = numeric(0))
  rownames(st) <- NULL
  list(keep = keep, markers = mk, dp = dp, r2 = r2, lod = lod, st = st, maf = maf, pos = pos, chr = ch,
       positions_used = used, max_dist_kb = if (is.finite(md)) md else 0, n_pairs = length(ut))
}

# ---- block definitions (Haploview 4.2 FindBlocks) ---------------------------
# st must contain every pair inside the maximum distance (hv_ld_engine()$st).
.hv_lookup <- function(st, n) { L <- matrix(NA_integer_, n, n); L[cbind(st$i, st$j)] <- seq_len(nrow(st)); L }
.hv_same <- function(ch, n) if (is.null(ch)) rep("1", n) else ch

#' Gabriel et al. (2002) as in Haploview: confidence-interval blocks.
#' A pair is in "strong LD" when its D' interval is [>= low_ci, >= 0.98] and
#' shows "historical recombination" when the upper bound is < 0.90; a region is
#' a block when at least 95% of its informative pairs are strong. Markers with
#' MAF < 0.05 are not used.
hv_blocks_gabriel <- function(st, pos, maf, chr = NULL, cutLowCI = .70, cutHighCI = .98, mafThresh = .05,
                              cutLowCIVar = c(0, 0, .8, .5, .5), maxDist = c(0, 0, 20000, 30000, 1e6),
                              recHighCI = .9, informFrac = .95) {
  n <- length(pos); L <- .hv_lookup(st, n); skip <- !(maf >= mafThresh); ch <- .hv_same(chr, n)
  cand <- which(!skip[st$i] & !skip[st$j] & st$LOD >= -90 & st$CIhi >= cutHighCI & st$CIlow >= cutLowCI)
  if (!length(cand)) return(list())
  cand <- data.frame(i = st$i[cand], j = st$j[cand], sep = abs(pos[st$j[cand]] - pos[st$i[cand]]))
  cand <- cand[order(cand$sep), , drop = FALSE]; cand <- cand[nrow(cand):1, , drop = FALSE]  # largest first
  inb <- logical(n); blocks <- list()
  for (r in seq_len(nrow(cand))) {
    a <- cand$i[r]; b <- cand$j[r]; sep <- cand$sep[r]
    if (inb[a] || inb[b]) next
    mk <- (a:b)[!skip[a:b]]; nm <- length(mk)
    if (nm < 4 && sep > maxDist[nm + 1]) next
    strong <- 0; rec <- 0
    for (x in mk[-1]) for (y in mk[mk < x]) {
      kx <- L[y, x]; if (is.na(kx)) next
      lodv <- st$LOD[kx]; lo <- st$CIlow[kx]; hi <- st$CIhi[kx]
      if (lodv < -90) next
      if (lodv == 0 && lo == 0 && hi == 0) next
      if (nm < 5) { if (lo > cutLowCIVar[nm + 1] && hi >= cutHighCI) strong <- strong + 1 }
      else if (lo > cutLowCI && hi >= cutHighCI) strong <- strong + 1
      if (hi < recHighCI) rec <- rec + 1
    }
    if (nm > 3) { if (strong + rec < 6) next } else if (nm > 2) { if (strong + rec < 3) next } else if (strong + rec < 1) next
    if (strong / (strong + rec) > informFrac) { blocks[[length(blocks) + 1L]] <- mk; inb[a:b] <- TRUE }
  }
  if (!length(blocks)) return(list())
  blocks[order(vapply(blocks, `[`, 1, 1))]
}

#' Solid spine of LD as in Haploview: the first and last marker of a block are
#' in strong LD (D' >= spine_dp) with every marker in between (one weak pair
#' along the spine is tolerated, as in Haploview).
hv_blocks_spine <- function(st, pos, chr = NULL, max_dist_kb = 500, spineDP = .8) {
  n <- length(pos); L <- .hv_lookup(st, n); ch <- .hv_same(chr, n)
  md <- suppressWarnings(as.numeric(max_dist_kb)); if (!is.finite(md) || md <= 0) md <- Inf
  stat <- function(a, b) { kx <- L[a, b]; if (is.na(kx)) NA_real_ else st$Dp[kx] }
  rowlen <- function(a) {
    if (a >= n) return(1)
    w <- which(ch[(a + 1):n] == ch[a] & abs(pos[(a + 1):n] - pos[a]) <= md * 1000)
    1 + if (length(w)) max(w) else 0
  }
  blocks <- list(); i <- 1L; last <- 0L
  while (i <= n) {
    mis <- 0; end <- 0L; j <- i + 1L; len <- rowlen(i)
    while (j < i + len) {
      d <- stat(i, j)
      if (!is.na(d)) {
        if (d < spineDP) { if (mis < 1) { mis <- mis + 1; end <- j } else { end <- j - 1L; break } } else end <- j
      }
      j <- j + 1L
    }
    k <- end
    while (k > i) {
      for (m in i:(k - 1L)) {
        d <- stat(m, k); if (is.na(d)) next
        if (d < spineDP) { if (mis >= 1) break; mis <- mis + 1 }
        last <- m + 1L
      }
      if (last == k) { blocks[[length(blocks) + 1L]] <- i:k; i <- k; break }
      k <- k - 1L
    }
    i <- i + 1L
  }
  blocks
}

#' Four-gamete rule (Wang et al. 2002) as in Haploview: a block has no pair of
#' markers in which all four two-marker haplotypes reach the cutoff (1%).
hv_blocks_4gamete <- function(st, pos, chr = NULL, cutoff = .01) {
  n <- length(pos); L <- .hv_lookup(st, n)
  thr <- cutoff + 1e-8
  four <- ((st$f00 > thr) + (st$f01 > thr) + (st$f10 > thr) + (st$f11 > thr)) > 3
  four[is.na(four)] <- FALSE
  cand <- which(!four)
  if (!length(cand)) return(list())
  cand <- data.frame(i = st$i[cand], j = st$j[cand], sep = st$j[cand] - st$i[cand] - 1L)
  cand <- cand[order(cand$sep), , drop = FALSE]; cand <- cand[nrow(cand):1, , drop = FALSE]
  inb <- logical(n); blocks <- list()
  for (r in seq_len(nrow(cand))) {
    a <- cand$i[r]; b <- cand$j[r]; if (inb[a] || inb[b]) next
    good <- TRUE
    for (x in (a + 1L):b) {
      for (y in a:(x - 1L)) { kx <- L[y, x]; if (!is.na(kx) && four[kx]) { good <- FALSE; break } }
      if (!good) break
    }
    if (good) { blocks[[length(blocks) + 1L]] <- a:b; inb[a:b] <- TRUE }
  }
  if (!length(blocks)) return(list())
  blocks[order(vapply(blocks, `[`, 1, 1))]
}

# ---- Haploview standard colour scheme ----------------------------------------
#' Haploview "Standard (D'/LOD)": LOD >= 2 (strictly > 2) and D' = 1 bright red,
#' LOD > 2 and D' < 1 shades of pink/red, LOD <= 2 and D' = 1 blue,
#' LOD <= 2 and D' < 1 white. Vectorised over dp, lod.
hv_color_standard <- function(dp, lod) {
  out <- rep("#FFFFFF", length(dp))
  ok <- is.finite(dp) & is.finite(lod)
  hiL <- ok & lod > 2
  pink <- hiL & dp < .5
  out[pink] <- grDevices::rgb(255, 224, 224, maxColorValue = 255)
  red <- hiL & !(dp < .5)
  if (any(red)) {
    v <- pmax(0, pmin(255, as.integer(446 * (1 - dp[red]))))
    out[red] <- grDevices::rgb(255, v, v, maxColorValue = 255)
  }
  blue <- ok & !(lod > 2) & dp > .99
  out[blue] <- grDevices::rgb(192, 192, 240, maxColorValue = 255)
  out
}

# Rough peak working memory (bytes) of hv_ld_engine() plus the figure for
# n individuals x k markers (used to refuse a job that would not fit).
hv_estimate_memory_bytes <- function(n, k) {
  n <- as.numeric(n); k <- as.numeric(k)
  140 * k^2 + 40 * n * k + 150e6
}
