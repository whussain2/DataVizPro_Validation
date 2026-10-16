# Install the R packages used by the scripts (R >= 4.5)
cran <- c("data.table", "mgcv", "spdep", "sp", "gstat", "SpATS", "clue", "FactoMineR", "irlba",
          "rrBLUP", "cluster", "MASS", "ggplot2", "patchwork")
install.packages(setdiff(cran, rownames(installed.packages())))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("LEA", "snpStats"), update = FALSE, ask = FALSE)
