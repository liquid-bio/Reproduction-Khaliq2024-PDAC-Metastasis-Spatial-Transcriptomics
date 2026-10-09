# 01 - Load Visium sections, spot-level QC, standard scRNA-style preprocessing
# Run from the repository root: Rscript R/01_qc.R
library(Seurat)
library(Matrix)
library(hdf5r)
library(jsonlite)
library(dplyr)

out <- "results/01_qc"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
samples <- read.csv("data/samples.csv")

# X in the h5ad is a spots x genes CSR matrix; reading it as CSC gives genes x spots
read_h5ad_counts <- function(f) {
  h <- H5File$new(f, "r")
  m <- sparseMatrix(i = h[["X/indices"]][] + 1, p = h[["X/indptr"]][], x = h[["X/data"]][],
                    dims = c(length(h[["var/_index"]][]), length(h[["obs/_index"]][])),
                    dimnames = list(h[["var/_index"]][], h[["obs/_index"]][]))
  h$close_all()
  m
}

# tissue_positions files come with a header in some sections and without one (*_list.csv) in others
read_positions <- function(dir) {
  f <- list.files(dir, "tissue_positions", full.names = TRUE)
  p <- read.csv(f, header = grepl("^barcode", readLines(f, 1)))
  names(p) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
  p
}

counts <- list(); spots <- list()
for (s in samples$sample) {
  d <- file.path("data/visium", s)
  m <- read_h5ad_counts(file.path(d, "expression.h5ad"))
  pos <- read_positions(file.path(d, "spatial")) %>% filter(in_tissue == 1)
  sf <- fromJSON(list.files(file.path(d, "spatial"), "scalefactors", full.names = TRUE))
  # Visium hexagonal grid: 6 neighbours at (row, col +/- 2) and (row +/- 1, col +/- 1).
  # Spots with < 6 in-tissue neighbours sit on the tissue edge.
  key <- paste(pos$array_row, pos$array_col)
  nb <- rowSums(sapply(list(c(0, 2), c(0, -2), c(1, 1), c(1, -1), c(-1, 1), c(-1, -1)),
                       function(o) paste(pos$array_row + o[1], pos$array_col + o[2]) %in% key))
  spots[[s]] <- pos %>% mutate(sample = s, n_neighbors = nb,
                               kept_by_authors = barcode %in% colnames(m),
                               x = pxl_col * sf$tissue_hires_scalef, y = pxl_row * sf$tissue_hires_scalef)
  colnames(m) <- paste(s, colnames(m), sep = "_")
  counts[[s]] <- m
}
spots <- bind_rows(spots) %>% left_join(samples, by = "sample") %>%
  mutate(spot_id = paste(sample, barcode, sep = "_"))

# Spots before/after the authors' QC (compare with Supplementary Table 2)
qc_count <- spots %>% group_by(sample, short, group) %>%
  summarise(in_tissue = n(), kept = sum(kept_by_authors), removed = in_tissue - kept, .groups = "drop")
write.csv(qc_count, file.path(out, "spots_per_sample.csv"), row.names = FALSE)

obj <- CreateSeuratObject(do.call(cbind, counts), project = "PDAC_Visium")
meta <- spots %>% filter(kept_by_authors) %>% select(-barcode, -in_tissue, -kept_by_authors)
obj <- AddMetaData(obj, data.frame(meta[match(colnames(obj), meta$spot_id), ], row.names = colnames(obj)))

# The FFPE probe panel contains no MT- or RPS/RPL genes, so only haemoglobin can be computed
cat("MT genes:", sum(grepl("^MT-", rownames(obj))), " RP genes:", sum(grepl("^RP[SL]", rownames(obj))), "\n")
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[AB]")
# Paper criterion: < 200 genes or UMIs = low quality. Flag only, do not remove.
obj$low_quality <- obj$nFeature_RNA < 200 | obj$nCount_RNA < 200

# Edge fraction of removed vs kept spots
spots <- spots %>% left_join(obj@meta.data %>% transmute(spot_id, nCount_RNA, nFeature_RNA, percent.hb, low_quality), by = "spot_id")
write.csv(spots, file.path(out, "spots_all_in_tissue.csv"), row.names = FALSE)
edge_tab <- spots %>% mutate(status = ifelse(kept_by_authors, ifelse(low_quality, "kept_low", "kept_ok"), "removed")) %>%
  group_by(status) %>% summarise(n = n(), frac_edge = mean(n_neighbors < 6), median_neighbors = median(n_neighbors))
write.csv(edge_tab, file.path(out, "edge_by_status.csv"), row.names = FALSE)
print(edge_tab)

# Standard scRNA-seq workflow, used to show that spots group by section rather than by tissue type
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)
obj <- RunUMAP(obj, dims = 1:30)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.3)

saveRDS(obj, file.path(out, "visium_qc.rds"))
