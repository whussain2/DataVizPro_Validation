# ==============================================================================
# 09_phenotypes.R - Case study 3 (Supplemental Figure 5) and the PCA / regression rows of Supplemental Table 4
# 18 traits of the rice diversity panel (Zhao et al. 2011)
# ==============================================================================
source("scripts/00_setup.R")
traits <- c("Flowering.time.at.Arkansas", "Plant.height", "Flag.leaf.length", "Flag.leaf.width",
            "Panicle.number.per.plant", "Panicle.length", "Primary.panicle.branch.number",
            "Seed.number.per.panicle", "Florets.per.panicle", "Panicle.fertility", "Seed.length",
            "Seed.width", "Seed.volume", "Seed.surface.area", "Brown.rice.seed.length",
            "Amylose.content", "Alkali.spreading.value", "Protein.content")
all <- read.csv(file.path(DATA_DIR, "rice", "rice_phenotypes.csv"), row.names = 1)
ph <- all[, traits]
cc <- ph[complete.cases(ph), ]

# ---- Correlations (DataVizPro compute_correlation, Benjamini-Hochberg)
co <- compute_correlation(ph, method = "pearson", p_adjust = "BH")
ut <- upper.tri(co$correlation)
pairs <- data.frame(trait_a = rownames(co$correlation)[row(co$correlation)[ut]],
                    trait_b = colnames(co$correlation)[col(co$correlation)[ut]],
                    r = co$correlation[ut], p_BH = co$p_value[ut])
pairs <- pairs[order(-abs(pairs$r)), ]
set.seed(1)
boot_ci <- function(a, b) {
  d <- ph[complete.cases(ph[, c(a, b)]), c(a, b)]
  quantile(replicate(1000, { i <- sample(nrow(d), replace = TRUE); cor(d[i, 1], d[i, 2]) }), c(.025, .975))
}
top <- head(pairs, 6)
ci <- t(mapply(boot_ci, top$trait_a, top$trait_b)); top$ci_low <- ci[, 1]; top$ci_high <- ci[, 2]
write.csv(pairs, file.path(TAB_DIR, "case3_trait_correlations.csv"), row.names = FALSE)
write.csv(top, file.path(TAB_DIR, "case3_top_correlations.csv"), row.names = FALSE)

# ---- PCA: DataVizPro (FactoMineR) versus prcomp
pca <- compute_pca(cc, n_components = 5, scale = TRUE, method = "factominer")
pr  <- prcomp(cc, center = TRUE, scale. = TRUE)
S_dvp <- as.matrix(pca$scores)[, 1:5]; S_ref <- pr$x[, 1:5]
for (k in 1:5) if (cor(S_dvp[, k], S_ref[, k]) < 0) S_dvp[, k] <- -S_dvp[, k]
a <- cmp(as.vector(S_dvp), as.vector(S_ref))
table4_row("PCA scores", "prcomp (base R)", "Rice phenotypes (359 x 18)", a["n"], a["r"], a["max_abs_diff"],
           "differences only from 1/n vs 1/(n-1) standardisation")
grp <- readRDS(file.path(TMP_DIR, "ancestry_groups.rds"))
scores <- data.frame(ID = rownames(cc), pr$x[, 1:2], group = grp[rownames(cc)])
write.csv(scores, file.path(TAB_DIR, "case3_phenotype_pcs.csv"), row.names = FALSE)
write.csv(round(pr$rotation[, 1:2], 3), file.path(TAB_DIR, "case3_pca_loadings.csv"))

# ---- Clustering: traits on 1 - |r| (average linkage); accessions Euclidean (Ward.D2)
library(cluster)
best_k <- function(d, h, ks) ks[which.max(sapply(ks, function(k) mean(silhouette(cutree(h, k), d)[, 3])))]
dt <- as.dist(1 - abs(cor(cc))); ht <- hclust(dt, "average"); kt <- best_k(dt, ht, 2:6)
da <- dist(scale(cc));           ha <- hclust(da, "ward.D2"); ka <- best_k(da, ha, 2:8)
trait_cl <- cutree(ht, kt)
write.csv(data.frame(trait = names(trait_cl), cluster = trait_cl), file.path(TAB_DIR, "case3_trait_clusters.csv"), row.names = FALSE)
acc_cl <- cutree(ha, ka)

summ <- data.frame(statistic = c("Accessions with complete records", "Trait pairs", "Pairs significant (BH 0.05)",
                                 "PC1 variance (%)", "PC2 variance (%)", "R2 of PC1 on ancestry group",
                                 "R2 of PC2 on ancestry group", "Trait clusters", "Trait cophenetic correlation",
                                 "Accession clusters", "Accession cophenetic correlation"),
                   value = c(nrow(cc), sum(ut), sum(pairs$p_BH < 0.05),
                             100 * (pr$sdev^2 / sum(pr$sdev^2))[1:2],
                             summary(lm(PC1 ~ group, scores))$r.squared, summary(lm(PC2 ~ group, scores))$r.squared,
                             kt, cor(cophenetic(ht), dt), ka, cor(cophenetic(ha), da)))
summ$value <- round(summ$value, 3)
write.csv(summ, file.path(TAB_DIR, "case3_summary.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(table(cluster = acc_cl, group = grp[names(acc_cl)])),
          file.path(TAB_DIR, "case3_accession_clusters_by_group.csv"))

# ---- Regression: DataVizPro fit_regression_model versus statsmodels (Python)
m1 <- fit_regression_model(all, "Plant.height", c("Flowering.time.at.Arkansas", "Panicle.length", "Seed.number.per.panicle"), "linear")
m2 <- fit_regression_model(all, "Awn.presence", c("Seed.length", "Plant.height", "Flag.leaf.width"), "logistic")
write.csv(m1$data[, c(m1$response, m1$predictors)], file.path(TMP_DIR, "ols.csv"), row.names = FALSE)
write.csv(m2$data[, c(m2$response, m2$predictors)], file.path(TMP_DIR, "logit.csv"), row.names = FALSE)
sm_file <- file.path(TMP_DIR, "statsmodels_coefs.csv")
if (nzchar(Sys.which("python3"))) system2("python3", c("python/references.py", "statsmodels"))
if (file.exists(sm_file)) {
  s <- read.csv(sm_file)
  a <- cmp(c(coef(m1$model), coef(m2$model)), c(s$ols, s$logit))
  table4_row("Regression coefficients", "statsmodels 0.15 (OLS, logistic)", "Rice phenotypes",
             a["n"], a["r"], a["max_abs_diff"], "linear and logistic models")
} else message("statsmodels output missing: run python/references.py statsmodels, then re-run this script")
print(summ)
