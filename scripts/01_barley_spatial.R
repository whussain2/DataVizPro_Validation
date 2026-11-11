# ==============================================================================
# 01_barley_spatial.R - Case study 1 (Figure 5) and the SpATS row of Supplemental Table 4
# Data: barley row-column trial of Durban et al. (2003), 272 genotypes x 2 reps,
#       544 plots on a 16 x 34 grid (agridat::durban.rowcol)
# ==============================================================================
source("scripts/00_setup.R")
suppressPackageStartupMessages(library(SpATS))

dat <- read.csv(file.path(DATA_DIR, "barley", "Durban2003Barley.csv"))

# ---- DataVizPro field-spatial module: genotype fixed, rep random, P-spline trend
field <- fs_build_field(dat, x_col = "bed", y_col = "row", value_col = "yield",
                        geno_col = "gen", block_cols = "rep")
ps <- fs_default_settings()
trend <- fs_fit_trend(field, ps)
m <- trend$model

spatial_share <- 100 * var(m$spatial) / var(m$data$z)
delta_aic     <- AIC(m$base_fit) - AIC(m$fit)

# Spatial diagnostics on residuals after design effects only, and after design + spatial
ps_design <- ps; ps_design$diag_scope <- "design"
ps_full   <- ps; ps_full$diag_scope   <- "full"
moran_design <- fs_morans_i(field, ps_design)
moran_full   <- fs_morans_i(field, ps_full)
lisa_design  <- fs_lisa(field, ps_design)$data
lisa_full    <- fs_lisa(field, ps_full)$data

# Semivariograms of the design-adjusted residuals (directional: 0 deg = rows, 90 deg = beds)
pts <- fs_diag_points(field, ps_design)
spdf <- data.frame(x = pts$x, y = pts$y, z = pts$z); sp::coordinates(spdf) <- ~ x + y
v_dir  <- gstat::variogram(z ~ 1, data = spdf, alpha = c(0, 90))
v_omni <- gstat::variogram(z ~ 1, data = spdf)
v_fit  <- gstat::fit.variogram(v_omni, gstat::vgm(var(pts$z) * 0.7, "Sph", max(v_omni$dist) / 2, var(pts$z) * 0.1))

# Genotype rankings: raw means versus spatially adjusted means (DataVizPro "adjusted" = y - spatial)
d <- m$data; d$adjusted <- d$z - m$spatial
raw_mean <- tapply(d$z, d$geno, mean)
adj_mean <- tapply(d$adjusted, d$geno, mean)
top27 <- function(v) names(sort(v, decreasing = TRUE))[1:27]   # top 10% of 272

# Genotype BLUEs from the DataVizPro model (spatial and rep terms excluded)
newd <- data.frame(geno = levels(d$geno), x = mean(d$x), y = mean(d$y))
for (b in m$block_names) newd[[b]] <- levels(d[[b]])[1]
smooth_labels <- sapply(m$fit$smooth, `[[`, "label")
blue_dvp <- predict(m$fit, newdata = newd, exclude = grep("^s\\(blk|^te\\(", smooth_labels, value = TRUE))
names(blue_dvp) <- newd$geno

# ---- Reference: SpATS (genotype fixed, PSANOVA surface, rep fixed)
dat$rep <- factor(dat$rep)
sm <- SpATS(response = "yield", spatial = ~ PSANOVA(bed, row, nseg = c(17, 8)),
            genotype = "gen", genotype.as.random = FALSE, fixed = ~ rep,
            data = dat, control = list(tolerance = 1e-03, monitoring = 0))
bs <- predict(sm, which = "gen")
blue_spats <- setNames(bs$predicted.values, bs$gen)
agree <- cmp(blue_dvp[names(blue_spats)], blue_spats)
table4_row("Spatially adjusted genotype values", "SpATS 1.0 (PSANOVA)", "Barley trial",
           agree["n"], agree["r"], agree["max_abs_diff"],
           sprintf("Spearman %.3f", cor(blue_dvp[names(blue_spats)], blue_spats, method = "spearman")))

# ---- Summary table for Section 3.5
cl <- function(x) sum(x$cluster != "Not significant")
summary_tab <- data.frame(
  statistic = c("Correlation of yield with row position", "Correlation of yield with bed position",
                "Spatial variance share (%)", "Deviance explained (%)", "Delta AIC (no spatial - spatial)",
                "P-spline basis (k x k)",
                "Moran's I, design residuals", "Moran's I z, design residuals", "Moran's I P, design residuals",
                "Moran's I, design + spatial residuals", "Moran's I P, design + spatial residuals",
                "LISA significant plots (FDR 0.05), design", "LISA significant plots (FDR 0.05), design + spatial",
                "Semivariogram nugget", "Semivariogram partial sill", "Semivariogram range (plots)",
                "Spearman raw vs adjusted genotype means", "Top-10% genotypes retained (of 27)"),
  value = c(cor(dat$yield, dat$row), cor(dat$yield, dat$bed),
            spatial_share, 100 * summary(m$fit)$dev.expl, delta_aic, m$k[1],
            moran_design$I, moran_design$z, moran_design$p, moran_full$I, moran_full$p,
            cl(lisa_design), cl(lisa_full),
            v_fit$psill[1], v_fit$psill[2], v_fit$range[2],
            cor(raw_mean, adj_mean, method = "spearman"), length(intersect(top27(raw_mean), top27(adj_mean)))))
summary_tab$value <- signif(summary_tab$value, 4)
write.csv(summary_tab, file.path(TAB_DIR, "case1_barley_summary.csv"), row.names = FALSE)
write.csv(data.frame(gen = names(blue_spats), raw_mean = raw_mean[names(blue_spats)],
                     adjusted_mean = adj_mean[names(blue_spats)], blue_datavizpro = blue_dvp[names(blue_spats)],
                     blue_spats = blue_spats),
          file.path(TAB_DIR, "case1_barley_genotypes.csv"), row.names = FALSE)
write.csv(as.data.frame(v_dir)[, c("dist", "gamma", "np", "dir.hor")], file.path(TAB_DIR, "case1_barley_variograms.csv"), row.names = FALSE)
saveRDS(list(trend = trend$data, lisa_design = lisa_design), file.path(TMP_DIR, "barley_maps.rds"))
print(summary_tab)
