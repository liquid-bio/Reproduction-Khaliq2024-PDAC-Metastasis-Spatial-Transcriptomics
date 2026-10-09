# 03 - MISTy spatial interaction analysis on RCTD cell-type proportions
# Run from the repository root: Rscript R/03_misty.R
library(mistyR)
library(future)
library(dplyr)

out <- "results/03_misty"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
plan(multisession, workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16")))

spots <- read.csv("results/01_qc/spots_all_in_tissue.csv") %>% filter(kept_by_authors)
samples <- read.csv("data/samples.csv")

for (s in samples$sample) {
  w <- read.csv(file.path("results/02_rctd", paste0("weights_", s, ".csv")), row.names = 1, check.names = FALSE)
  colnames(w) <- make.names(colnames(w))
  p <- spots[match(rownames(w), spots$spot_id), ]
  # Array coordinates -> planar coordinates in units of the spot centre-to-centre distance (100 um)
  pos <- data.frame(x = p$array_col / 2, y = p$array_row * sqrt(3) / 2)
  # As in the paper: intraview + juxtaview (radius 5 spots) + paraview (15 spots)
  create_initial_view(w) %>%
    add_juxtaview(pos, neighbor.thr = 5) %>%
    add_paraview(pos, l = 15) %>%
    run_misty(results.folder = file.path(out, "runs", s))
}

res <- collect_results(file.path(out, "runs", samples$sample))
imp <- res$importances %>%
  mutate(sample = basename(sample)) %>%
  left_join(samples %>% select(sample, group), by = "sample")
write.csv(imp, file.path(out, "importances_per_sample.csv"), row.names = FALSE)
write.csv(res$improvements %>% mutate(sample = basename(sample)), file.path(out, "improvements_per_sample.csv"), row.names = FALSE)

# Median standardised importance across sections (all sections and per tissue group)
med <- imp %>% group_by(view, Predictor, Target) %>% summarise(median_importance = median(Importance, na.rm = TRUE), .groups = "drop")
med_grp <- imp %>% group_by(group, view, Predictor, Target) %>% summarise(median_importance = median(Importance, na.rm = TRUE), .groups = "drop")
write.csv(med, file.path(out, "importance_median_all.csv"), row.names = FALSE)
write.csv(med_grp, file.path(out, "importance_median_by_group.csv"), row.names = FALSE)
