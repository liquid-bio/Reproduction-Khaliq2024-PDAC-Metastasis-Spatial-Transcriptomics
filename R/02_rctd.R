# 02 - RCTD deconvolution (full mode) against the authors' integrated scRNA-seq reference
# (column celltype_new1 = Celltype_Tier2, 15 cell types)
# Run from the repository root: Rscript R/02_rctd.R [sample ...]   (no argument = all sections)
library(Seurat)
library(spacexr)
library(Matrix)

out <- "results/02_rctd"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(TRUE)
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))

sc <- readRDS("data/ref/scRNA-seq_Data_post_qc.rds")
sc_counts <- GetAssayData(sc, assay = "RNA", layer = "counts")
ref <- Reference(sc_counts, factor(setNames(sc$celltype_new1, colnames(sc))), colSums(sc_counts))
rm(sc, sc_counts); gc()

obj <- readRDS("results/01_qc/visium_qc.rds")
st_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
run_samples <- if (length(args)) args else unique(obj$sample)

for (s in run_samples) {
  cells <- colnames(obj)[obj$sample == s]
  puck <- SpatialRNA(data.frame(x = obj$x[cells], y = obj$y[cells], row.names = cells),
                     st_counts[, cells], colSums(st_counts[, cells]))
  rctd <- create.RCTD(puck, ref, max_cores = cores, UMI_min = 20)
  rctd <- run.RCTD(rctd, doublet_mode = "full")
  w <- as.matrix(normalize_weights(rctd@results$weights))
  write.csv(w, file.path(out, paste0("weights_", s, ".csv")))
}
