# Load the DataVizPro computation engine (extract of v1.1.1; haploview.R from v1.2) into the global environment.
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a
engine_dir <- if (exists("ENGINE_DIR")) ENGINE_DIR else "dvp_engine"
for (f in c("fieldspatial", "ps_io", "ps_inference", "ld_decay", "ld_heatmap",
            "genomics_heatmap", "pca_helpers", "pca", "regression", "correlation", "haploview")) {
  sys.source(file.path(engine_dir, paste0(f, ".R")), envir = globalenv(), keep.source = FALSE)
}
rm(f)
