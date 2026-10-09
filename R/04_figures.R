# 04 - Figures and summary tables (reads only the outputs of 01-03)
# Run from the repository root: Rscript R/04_figures.R
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(png)

fig <- "figures"
tab <- "results/tables"
dir.create(fig, showWarnings = FALSE)
dir.create(tab, recursive = TRUE, showWarnings = FALSE)
samples <- read.csv("data/samples.csv")
grp_col <- c(PT = "#d95f02", HM = "#7570b3", LNM = "#1b9e77", NP = "#e7298a")
theme_set(theme_classic(base_size = 9))

# Plot spots on the hires H&E image; df needs x, y in hires pixels
he <- function(s) readPNG(list.files(file.path("data/visium", s, "spatial"), "hires_image", full.names = TRUE))
spatial_plot <- function(df, s, value, title = value, size = 0.5, discrete = FALSE) {
  img <- he(s)
  d <- df[df$sample == s, ]
  pad <- 40
  p <- ggplot(d, aes(x, y)) +
    annotation_raster(img, xmin = 0, xmax = ncol(img), ymin = -nrow(img), ymax = 0) +
    geom_point(aes(colour = .data[[value]], y = -y), size = size, stroke = 0) +
    coord_equal(xlim = range(d$x) + c(-pad, pad), ylim = -rev(range(d$y)) + c(-pad, pad), expand = FALSE) +
    labs(title = title, colour = NULL) + theme_void(base_size = 8) +
    theme(legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))
  if (discrete) p else p + scale_colour_viridis_c(option = "magma")
}
he_only <- function(s, title) {
  img <- he(s)
  d <- spots[spots$sample == s, ]
  pad <- 40
  ggplot() + annotation_raster(img, xmin = 0, xmax = ncol(img), ymin = -nrow(img), ymax = 0) +
    geom_blank(aes(x, -y), d) +
    coord_equal(xlim = range(d$x) + c(-pad, pad), ylim = -rev(range(d$y)) + c(-pad, pad), expand = FALSE) +
    labs(title = title) + theme_void(base_size = 8) + theme(plot.title = element_text(hjust = 0.5))
}
short <- setNames(samples$short, samples$sample)

# ---------- 1. QC ----------
spots <- read.csv("results/01_qc/spots_all_in_tissue.csv")
kept <- spots %>% filter(kept_by_authors) %>%
  mutate(short = factor(short, samples$short[order(factor(samples$group, names(grp_col)), samples$short)]))

p1a <- kept %>% pivot_longer(c(nCount_RNA, nFeature_RNA)) %>%
  ggplot(aes(short, value, fill = group)) + geom_violin(scale = "width", linewidth = 0.2) +
  geom_hline(yintercept = 200, linetype = 2, colour = "red") +
  scale_y_log10() + scale_fill_manual(values = grp_col) + facet_wrap(~name, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(file.path(fig, "fig1a_qc_violin.png"), p1a, width = 7, height = 4, dpi = 300, bg = "white")

# Where the authors removed spots: the 4 sections with the largest removed fraction
cnt <- read.csv("results/01_qc/spots_per_sample.csv") %>% mutate(frac_removed = removed / in_tissue) %>% arrange(-frac_removed)
write.csv(cnt, file.path(tab, "spots_removed_per_section.csv"), row.names = FALSE)
top4 <- head(cnt$sample, 4)
spots$status <- ifelse(spots$kept_by_authors, "kept", "removed")
p1b <- wrap_plots(c(lapply(top4, function(s) he_only(s, short[s])),
                    lapply(top4, function(s) spatial_plot(spots, s, "status", "", discrete = TRUE) +
                             scale_colour_manual(values = c(kept = "grey70", removed = "red"))),
                    lapply(top4, function(s) spatial_plot(mutate(kept, log10_UMI = log10(nCount_RNA)), s, "log10_UMI", "") + scale_colour_viridis_c(option = "magma", limits = c(2.3, 4.5), oob = scales::squish))),
                  nrow = 3, guides = "collect")
ggsave(file.path(fig, "fig1b_removed_spots.png"), p1b, width = 8, height = 6.5, dpi = 300, bg = "white")

# Fraction of edge spots: removed vs kept
p1c <- spots %>% mutate(edge = n_neighbors < 6) %>% group_by(sample, status) %>% summarise(frac_edge = mean(edge), .groups = "drop") %>%
  ggplot(aes(status, frac_edge)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, size = 0.8) +
  labs(x = NULL, y = "fraction of edge spots (<6 neighbors)")
# UMI of kept spots on the edge vs in the interior
p1d <- kept %>% mutate(edge = ifelse(n_neighbors < 6, "edge", "interior")) %>%
  ggplot(aes(edge, nCount_RNA)) + geom_boxplot(outlier.size = 0.1) + scale_y_log10() + labs(x = NULL, y = "UMI per spot (kept)")
ggsave(file.path(fig, "fig1c_edge.png"), p1c + p1d, width = 4.5, height = 2.5, dpi = 300, bg = "white")

# Proposed per-section criterion: log10 UMI below median - 3 MAD within the section (flag only)
mad_flag <- kept %>% group_by(sample, short, group) %>%
  mutate(log_umi = log10(nCount_RNA), mad_outlier = log_umi < median(log_umi) - 3 * mad(log_umi)) %>% ungroup()
write.csv(mad_flag %>% group_by(short, group) %>%
            summarise(n = n(), umi_threshold = round(10^(median(log_umi) - 3 * mad(log_umi))), n_flagged = sum(mad_outlier),
                      frac_edge_flagged = mean(n_neighbors[mad_outlier] < 6), .groups = "drop") %>% arrange(-n_flagged),
          file.path(tab, "qc_mad_outliers_per_section.csv"), row.names = FALSE)
mad_flag$flag <- ifelse(mad_flag$mad_outlier, "3-MAD outlier", "pass")
p1f <- wrap_plots(lapply(c("HM3343", "HM3357"), function(s) spatial_plot(mad_flag, s, "flag", short[s], size = 0.6, discrete = TRUE) +
                           scale_colour_manual(values = c(pass = "grey70", "3-MAD outlier" = "red"))), guides = "collect")
ggsave(file.path(fig, "fig1d_mad_outliers.png"), p1f, width = 5, height = 2.4, dpi = 300, bg = "white")

# UMAP from the scRNA-style workflow: section effect
obj <- readRDS("results/01_qc/visium_qc.rds")
obj$short <- unname(short[obj$sample])
p1e <- DimPlot(obj, group.by = "group", cols = grp_col, raster = TRUE) + ggtitle("tissue group") |
  DimPlot(obj, group.by = "short", raster = TRUE) + ggtitle("section") + theme(legend.text = element_text(size = 5)) + guides(colour = guide_legend(ncol = 2, override.aes = list(size = 1)))
ggsave(file.path(fig, "fig1e_umap.png"), p1e, width = 8, height = 3.5, dpi = 300, bg = "white")

# ---------- 2. RCTD deconvolution ----------
w <- do.call(rbind, lapply(samples$sample, function(s)
  read.csv(file.path("results/02_rctd", paste0("weights_", s, ".csv")), row.names = 1, check.names = FALSE)))
ct <- colnames(w)
obj <- subset(obj, cells = rownames(w))
obj <- AddMetaData(obj, w)
obj$Fibroblasts <- obj$myCAF + obj$iCAF
write.csv(w, gzfile(file.path(tab, "rctd_weights.csv.gz")))

# Mean cell-type proportion per section (cf. paper ED Fig. 1c)
p2a <- obj@meta.data %>% group_by(short, group) %>% summarise(across(all_of(ct), mean), .groups = "drop") %>%
  pivot_longer(all_of(ct), names_to = "celltype", values_to = "prop") %>%
  mutate(short = factor(short, levels(kept$short))) %>%
  ggplot(aes(short, prop, fill = celltype)) + geom_col(width = 0.9) +
  facet_grid(~factor(group, names(grp_col)), scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(12, "Paired"), "grey40", "black", "gold")) +
  labs(x = NULL, y = "mean RCTD proportion", fill = NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.key.size = unit(0.3, "cm"))
ggsave(file.path(fig, "fig2a_rctd_bar.png"), p2a, width = 7, height = 3.2, dpi = 300, bg = "white")
write.csv(obj@meta.data %>% group_by(group) %>% summarise(across(all_of(ct), mean)), file.path(tab, "rctd_mean_by_group.csv"), row.names = FALSE)

# Spot-level correlation between marker genes and RCTD proportions (sanity check of the deconvolution)
markers <- c("CEACAM5", "S100P", "PRSS1", "CPA1", "ALB", "APOA1", "COL11A1", "POSTN", "CFD", "C3",
             "RGS5", "NOTCH3", "PECAM1", "VWF", "MS4A1", "CD79A", "IL7R", "CD8A", "NKG7", "MKI67",
             "CD1C", "LAMP3", "SPP1", "C1QA", "C1QB", "FCN1", "S100A8")
expr <- t(as.matrix(GetAssayData(obj, layer = "data")[markers, ]))
ct_order <- c("Tumor Epithelial cells", "Normal Epithelial cells", "Hepatocytes", "myCAF", "iCAF", "PVL", "Endothelial cells",
              "B cells", "CD4+ cells", "CD8-NK cells", "Proliferative T cells", "DCs", "SPP1-TAM", "C1Q-TAM", "FCN1-TAM")
cm <- cor(expr, as.matrix(obj@meta.data[, ct]), method = "spearman")
p2b <- as.data.frame(as.table(cm)) %>%
  mutate(Var1 = factor(Var1, rev(markers)), Var2 = factor(Var2, ct_order)) %>%
  ggplot(aes(Var2, Var1, fill = Freq)) + geom_tile() +
  scale_fill_gradient2(low = "#2166ac", high = "#b2182b", name = "Spearman") +
  labs(x = "RCTD cell type", y = "marker gene") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig, "fig2b_marker_cor.png"), p2b, width = 5, height = 4.8, dpi = 300, bg = "white")

# Counterpart of paper Fig. 1c-e: the same 4 representative sections
rep4 <- c(PT2 = "HM3334", HM5 = "HM3341", LNM6 = "HM3343", NP2 = "HM3333")
# Fges gene sets from Bagaev et al. 2021 (CAF / tumour proliferation rate / B cells)
fges <- list(CAF = c("COL1A1", "COL1A2", "COL5A1", "COL6A1", "COL6A2", "COL6A3", "COL11A1", "LUM", "ACTA2", "FAP",
                     "PDGFRA", "PDGFRB", "FGF2", "CXCL12", "LRP1", "MMP2", "CD248", "FBLN1", "ELN", "MFAP5"),
             Tumor_prolif = c("MKI67", "ESCO2", "CETN3", "CDK2", "CCND1", "CCNE1", "AURKA", "AURKB", "E2F1", "MYBL2",
                              "BUB1", "PLK1", "CCNB1", "MCM2", "MCM6"),
             B_cells = c("CD19", "MS4A1", "TNFRSF13C", "CD22", "CR2", "FCRL5", "CD79A", "CD79B", "BLK", "PAX5",
                         "TNFRSF17", "STAP1"))
obj <- AddModuleScore(obj, fges, name = "fges")
names(obj@meta.data)[grep("^fges", names(obj@meta.data))] <- paste0(names(fges), "_score")
md <- cbind(obj@meta.data, t(as.matrix(GetAssayData(obj, layer = "data")[c("COL1A1", "KRT8", "KRT18", "EPCAM", "CD52", "CD79A"), ])))
epi_gene <- c(PT2 = "KRT18", HM5 = "KRT8", LNM6 = "KRT8", NP2 = "EPCAM")
imm_gene <- c(PT2 = "CD52", HM5 = "CD52", LNM6 = "CD79A", NP2 = "CD79A")
row_plots <- function(nm) {
  s <- rep4[[nm]]
  list(he_only(s, nm), spatial_plot(md, s, "COL1A1"), spatial_plot(md, s, epi_gene[[nm]]), spatial_plot(md, s, imm_gene[[nm]]),
       spatial_plot(md, s, "Fibroblasts"), spatial_plot(md, s, "Tumor Epithelial cells", "Tumor epi."), spatial_plot(md, s, "B cells"),
       spatial_plot(md, s, "CAF_score", "Fges CAF"), spatial_plot(md, s, "Tumor_prolif_score", "Fges prolif."), spatial_plot(md, s, "B_cells_score", "Fges B"))
}
p2c <- wrap_plots(do.call(c, lapply(names(rep4), row_plots)), ncol = 10) & theme(legend.position = "none", plot.title = element_text(size = 6))
ggsave(file.path(fig, "fig2c_insitu.png"), p2c, width = 10, height = 4.6, dpi = 300, bg = "white")

# ---------- 3. MISTy ----------
med <- read.csv("results/03_misty/importance_median_all.csv")
impr <- read.csv("results/03_misty/improvements_per_sample.csv")
lab <- function(x) unname(setNames(ct_order, make.names(ct_order))[x])   # undo make.names() on cell-type names
# Contribution of each view: intraview R2 and the gain from juxta/para (cf. paper ED Fig. 2a)
p3a <- impr %>% filter(measure %in% c("intra.R2", "gain.R2")) %>% group_by(target, measure) %>%
  summarise(value = median(value), .groups = "drop") %>%
  ggplot(aes(lab(target), value, fill = measure)) + geom_col() + coord_flip() +
  scale_fill_manual(values = c(intra.R2 = "grey60", gain.R2 = "#d95f02"), labels = c(gain.R2 = "gain from juxta+para", intra.R2 = "intraview R2")) +
  labs(x = NULL, y = "median R2 (%) across sections", fill = NULL) + theme(legend.position = "bottom")
heat <- function(v) med %>% filter(view == v) %>%
  ggplot(aes(lab(Predictor), lab(Target), fill = pmax(median_importance, 0))) + geom_tile() +
  scale_fill_gradient(low = "white", high = "#b2182b", limits = c(0, max(med$median_importance, na.rm = TRUE)), name = "median\nimportance") +
  labs(title = v, x = "predictor", y = "target") + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 6))
p3b <- heat("intra") + heat("juxta.5") + heat("para.15") + plot_layout(guides = "collect")
ggsave(file.path(fig, "fig3a_view_contribution.png"), p3a, width = 3.5, height = 3.5, dpi = 300, bg = "white")
ggsave(file.path(fig, "fig3b_importance_heatmaps.png"), p3b, width = 12, height = 4.3, dpi = 300, bg = "white")

# Intraview interaction network (cf. paper Fig. 2a): edges with median importance > 1, circular layout
edges <- med %>% filter(view == "intra", median_importance > 1)
nodes <- data.frame(name = sort(unique(c(edges$Predictor, edges$Target))))
nodes$a <- seq(0, 2 * pi, length.out = nrow(nodes) + 1)[-1]
nodes$x <- cos(nodes$a); nodes$y <- sin(nodes$a)
e <- edges %>% left_join(nodes, by = c(Predictor = "name")) %>% left_join(nodes, by = c(Target = "name"), suffix = c("", "end"))
p3c <- ggplot() + geom_curve(aes(x, y, xend = xend, yend = yend, linewidth = median_importance), e, colour = "grey50", alpha = 0.6, curvature = 0.15) +
  geom_point(aes(x, y), nodes, size = 4, colour = "#d95f02") + geom_text(aes(x * 1.1, y * 1.1, label = lab(name), hjust = ifelse(x < -0.01, 1, ifelse(x > 0.01, 0, 0.5))), nodes, size = 2.4) +
  scale_linewidth(range = c(0.3, 2), name = "median\nimportance") + coord_equal(xlim = c(-2.3, 2.6), ylim = c(-1.2, 1.2)) + theme_void(base_size = 8)
ggsave(file.path(fig, "fig3c_intra_network.png"), p3c, width = 4.5, height = 3.8, dpi = 300, bg = "white")

# Histology: distance of myCAF-/iCAF-rich spots to tumour (PT sections)
# Distance from each spot to the nearest tumour-rich spot (tumour epithelial > 0.5), in units of 100 um
pt <- obj@meta.data %>% filter(group == "PT") %>% mutate(px = array_col / 2, py = array_row * sqrt(3) / 2)
dist_tab <- pt %>% group_by(sample) %>% group_modify(function(d, k) {
  tum <- d[d$`Tumor Epithelial cells` > 0.5, ]
  dmin <- sapply(seq_len(nrow(d)), function(i) min(sqrt((tum$px - d$px[i])^2 + (tum$py - d$py[i])^2)))
  bind_rows(tibble(type = "myCAF-rich", dist = dmin[d$myCAF > 0.3]), tibble(type = "iCAF-rich", dist = dmin[d$iCAF > 0.3]))
}) %>% ungroup()
write.csv(dist_tab %>% group_by(sample, type) %>% summarise(n = n(), median_dist = median(dist), .groups = "drop"),
          file.path(tab, "caf_tumor_distance.csv"), row.names = FALSE)
p3d <- dist_tab %>% group_by(sample, type) %>% summarise(n = n(), median_dist = median(dist), .groups = "drop") %>%
  group_by(sample) %>% filter(all(n >= 20), n() == 2) %>% ungroup() %>%   # keep sections with >= 20 spots of both types
  ggplot(aes(type, median_dist)) + geom_boxplot(outlier.shape = NA) + geom_line(aes(group = sample), colour = "grey70") + geom_point(size = 0.8) +
  labs(x = NULL, y = "median distance to tumor-rich spot\n(×100 µm), per PT section")
s <- "HM3360"  # PT12
p3e <- wrap_plots(he_only(s, "PT12 H&E"), spatial_plot(md, s, "Tumor Epithelial cells", "Tumor epi."), spatial_plot(md, s, "myCAF"), spatial_plot(md, s, "iCAF"), p3d,
                  nrow = 1, widths = c(1, 1, 1, 1, 0.9))
ggsave(file.path(fig, "fig3d_caf_histology.png"), p3e, width = 10, height = 2.6, dpi = 300, bg = "white")
