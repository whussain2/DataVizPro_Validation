# ==============================================================================
# 10_figures.R - repository figures (results/figures) from the outputs of 01-09
# Colours follow the DataVizPro / PlantNura brand palette
# ==============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

BRAND <- c(purple = "#7B3FE4", green = "#2E8B3E", amber = "#F59E0B", blue = "#2F7BD9", magenta = "#D63F8C")
INK <- "#1F2937"; MUTED <- "#6B7280"; GRID <- "#E5E7EB"
GROUP_COLS <- c(G1 = BRAND[["purple"]], G2 = BRAND[["green"]], G3 = BRAND[["amber"]],
                G4 = BRAND[["blue"]], G5 = BRAND[["magenta"]], admixed = "#B8BEC8")
theme_dvp <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(text = element_text(colour = INK), plot.title = element_text(face = "bold", size = base + 2),
          plot.subtitle = element_text(colour = MUTED), axis.text = element_text(colour = MUTED),
          panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = GRID, linewidth = 0.3),
          strip.text = element_text(face = "bold", colour = INK, hjust = 0),
          plot.background = element_rect(fill = "white", colour = NA), legend.position = "bottom")
}
save_fig <- function(p, name, w, h) ggsave(file.path(FIG_DIR, name), p, width = w, height = h, dpi = 200, bg = "white")
set.seed(7)

# ---- 1. Validation: DataVizPro versus independent references -----------------
t4 <- read.csv(file.path(TAB_DIR, "table4_agreement.csv"))
take <- function(d, n = 6000) d[sample(nrow(d), min(n, nrow(d))), ]
# Relationship matrix
K <- readRDS(file.path(TMP_DIR, "grm.rds")); prep <- readRDS(file.path(TMP_DIR, "prep_all.rds"))
A <- rrBLUP::A.mat(prep$geno - 1, min.MAF = 0, max.missing = 1, impute.method = "mean", shrink = FALSE)
ut <- which(upper.tri(K))
d1 <- take(data.frame(x = A[ut], y = K[ut], panel = "Genomic relationships\nvs rrBLUP"))
# Ancestry
em <- readRDS(file.path(TMP_DIR, "em_k5.rds"))$q_by_k[["5"]]$best
fm <- as.matrix(read.table(file.path(TMP_DIR, "fm_k5.K5.s1.Q")))
s <- cor(em, fm); fm <- fm[, as.integer(clue::solve_LSAP(s - min(s), maximum = TRUE))]
d2 <- data.frame(x = as.vector(fm), y = as.vector(em), panel = "Ancestry coefficients (EM)\nvs fastmixture")
# LD r2 on 500 chromosome-1 markers
map <- read_rice_map(); chr <- map$chr[match(colnames(prep$geno), map$snp)]
M <- prep$geno[, which(chr == 1)[1:500]]
r2_dvp <- ld_compute(M, "r2"); r2_cor <- suppressWarnings(cor(M, use = "pairwise.complete.obs")^2)
u2 <- which(upper.tri(r2_dvp))
d3 <- take(data.frame(x = r2_cor[u2], y = r2_dvp[u2], panel = "LD r²\nvs genotype correlation"))
# Field-spatial
b <- read.csv(file.path(TAB_DIR, "case1_barley_genotypes.csv"))
d4 <- data.frame(x = b$blue_spats, y = b$blue_datavizpro, panel = "Adjusted genotype values\nvs SpATS")
dv <- rbind(d1, d2, d3, d4)
lab <- aggregate(cbind(x, y) ~ panel, dv, function(v) min(v))
lab$r <- sapply(lab$panel, function(p) { z <- dv[dv$panel == p, ]; cor(z$x, z$y) })
lab$label <- sprintf("r = %.4f", lab$r)
p_val <- ggplot(dv, aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, colour = MUTED, linewidth = 0.4, linetype = "dashed") +
  geom_point(colour = BRAND[["purple"]], alpha = 0.35, size = 0.9, stroke = 0) +
  geom_label(data = lab, aes(label = label), x = -Inf, y = Inf, hjust = -0.08, vjust = 1.3,
             size = 3.3, linewidth = 0, fill = "white", colour = INK) +
  facet_wrap(~ panel, nrow = 1, scales = "free") +
  labs(title = "DataVizPro reproduces independent reference software",
       subtitle = "Each point is one value computed by the DataVizPro engine (y) and by the reference (x); dashed line = identity",
       x = "Reference implementation", y = "DataVizPro") + theme_dvp()
save_fig(p_val, "fig_validation.png", 13, 4.2)

# ---- 2. Barley field: raw, trend, adjusted ------------------------------------
bm <- readRDS(file.path(TMP_DIR, "barley_maps.rds"))$trend
lims <- range(c(bm$observed, bm$adjusted))
tile <- function(col, title, lim, pal) {
  ggplot(bm, aes(x, y, fill = .data[[col]])) + geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradientn(colours = pal, limits = lim, name = NULL) + coord_equal(expand = FALSE) +
    labs(title = title, x = "Bed (column)", y = "Row") + theme_dvp(10) +
    theme(panel.grid = element_blank(), legend.key.width = unit(1.6, "cm"), legend.key.height = unit(0.25, "cm"))
}
seq_pal <- c("#F3EEFD", "#C9B3F5", "#9B6DEB", "#7B3FE4", "#4B1FA3")
div_pal <- c("#2F7BD9", "#A9C8F0", "#F3F4F6", "#E9A9C9", "#D63F8C")
tl <- max(abs(bm$fitted))
p_bar <- tile("observed", "Raw yield", lims, seq_pal) + tile("fitted", "Fitted P-spline trend", c(-tl, tl), div_pal) +
  tile("adjusted", "Spatially adjusted yield", lims, seq_pal) +
  plot_annotation(title = "Case study 1 · barley row–column trial (Durban et al., 2003)",
                  subtitle = "The spatial term explains 38.2% of yield variance (ΔAIC = 408); residual Moran's I falls from 0.255 to 0.037",
                  theme = theme_dvp())
save_fig(p_bar, "fig_barley_spatial.png", 13, 3.7)

# ---- 3. Manhattan plot: seed length --------------------------------------------
g <- read.csv(file.path(TMP_DIR, "gwas_Seed.length.csv"))
g <- g[order(g$CHR, g$BP), ]; g$logP <- -log10(g$P)
off <- cumsum(c(0, tapply(g$BP, g$CHR, max)[-12])); names(off) <- 1:12
g$pos <- g$BP + off[as.character(g$CHR)]
mids <- tapply(g$pos, g$CHR, function(v) mean(range(v)))
thr <- -log10(0.05 / nrow(g))
gs3 <- g[g$SNP == "wd3000590", ]
p_man <- ggplot(g, aes(pos, logP, colour = factor(CHR %% 2))) +
  geom_hline(yintercept = thr, colour = BRAND[["magenta"]], linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 0.7, alpha = 0.8, stroke = 0) +
  scale_colour_manual(values = c("0" = BRAND[["purple"]], "1" = "#9AA3AF"), guide = "none") +
  geom_point(data = gs3, colour = BRAND[["green"]], size = 2.6) +
  annotate("label", x = gs3$pos, y = gs3$logP + 0.9, label = "GS3", fontface = "italic", size = 3.4,
           linewidth = 0, fill = "white", colour = INK) +
  scale_x_continuous(breaks = mids, labels = 1:12, expand = expansion(0.01)) +
  labs(title = "Case study 2 · seed length GWAS, rice 44K panel",
       subtitle = sprintf("rrBLUP mixed model, 377 accessions, %s SNPs; dashed line = Bonferroni threshold", format(nrow(g), big.mark = ",")),
       x = "Chromosome", y = expression(-log[10](italic(P)))) + theme_dvp() + theme(panel.grid.major.x = element_blank())
save_fig(p_man, "fig_manhattan_seed_length.png", 11, 4)

# ---- 4. LD decay -----------------------------------------------------------------
ldd <- read.csv(file.path(TAB_DIR, "case2_rice_ld_decay.csv"))
curves <- do.call(rbind, lapply(seq_len(nrow(ldd)), function(i) {
  bn <- readRDS(file.path(TMP_DIR, sprintf("ld_decay_binned_%s.rds", ldd$group[i])))
  n <- 2 * ldd$accessions[i]; rho <- ldd$rho_datavizpro[i]; d <- seq(1, 1e6, length.out = 400); C <- rho * d
  fit <- ((10 + C) / ((2 + C) * (11 + C))) * (1 + ((3 + C) * (12 + 12 * C + C^2)) / (n * (2 + C) * (11 + C)))
  list(pts = data.frame(group = ldd$group[i], dist = bn$dist, r2 = bn$mean_r2),
       fit = data.frame(group = ldd$group[i], dist = d, r2 = fit))
}))
pts <- do.call(rbind, curves[, 1]); fits <- do.call(rbind, curves[, 2])
ldd$label <- sprintf("%s: half-decay %.0f kb", ldd$group, ldd$half_decay_kb_datavizpro)
p_ld <- ggplot() +
  geom_point(data = pts, aes(dist / 1e3, r2, colour = group), size = 1, alpha = 0.45, stroke = 0) +
  geom_line(data = fits, aes(dist / 1e3, r2, colour = group), linewidth = 1) +
  scale_colour_manual(values = GROUP_COLS, labels = setNames(ldd$label, ldd$group), name = NULL) +
  labs(title = "LD decay by ancestry group", subtitle = "Binned mean r² and fitted Hill–Weir curves (DataVizPro LD decay module)",
       x = "Distance (kb)", y = expression(r^2)) + theme_dvp()
save_fig(p_ld, "fig_ld_decay.png", 7, 4.4)

# ---- 5. Ancestry bar plot --------------------------------------------------------
q <- read.csv(file.path(TAB_DIR, "case2_rice_ancestry_K5.csv"))
qcols <- grep("^Cluster", names(q), value = TRUE)
q$order <- order(order(q$group, -apply(q[, qcols], 1, max)))
long <- do.call(rbind, lapply(seq_along(qcols), function(k)
  data.frame(acc = q$order, cluster = paste0("G", k), q = q[[qcols[k]]])))
p_q <- ggplot(long, aes(acc, q, fill = cluster)) + geom_col(width = 1) +
  scale_fill_manual(values = GROUP_COLS, name = "Ancestral cluster") +
  scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Ancestry proportions at K = 5", subtitle = "EM admixture engine; 413 accessions ordered by ancestry group",
       x = NULL, y = "Ancestry proportion") +
  theme_dvp() + theme(axis.text.x = element_blank(), panel.grid = element_blank())
save_fig(p_q, "fig_ancestry_K5.png", 11, 3.6)
message("Figures written to ", FIG_DIR)
