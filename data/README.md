# Input data

Raw data are not tracked in this repository. The scripts expect the following layout:

```
data/
├── samples.csv                      # tracked: section ID, GEO accession, tissue, paper section name, group
├── visium/<section>/expression.h5ad # spot x gene counts (one folder per section, 29 sections)
├── visium/<section>/spatial/        # *_tissue_positions*.csv, *_scalefactors_json.json, *_tissue_hires_image.png
└── ref/scRNA-seq_Data_post_qc.rds   # scRNA-seq reference (Seurat object)
```

Sources:

- Visium FFPE data from Khaliq et al., *Nat Genet* 2024, GEO [GSE272362](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272362). GEO only hosts images and spot positions; the per-section `expression.h5ad` count matrices were distributed with the course assignment. They contain the spots that passed the authors' QC (spot numbers match Supplementary Table 2 of the paper).
- `HM3335` (IU_PDA_HM2_2) is excluded: its count matrix in the course distribution is identical to `HM3336` (IU_PDA_HM2) and does not match its own tissue positions.
- scRNA-seq reference: Zenodo record [22695225](https://zenodo.org/records/22695225), file `scRNA-seq_Data_post_qc.rds` (1.18 GB, md5 `85b48bbceb1f9fd854fa5afb49e8a9d3`). Cell-type labels are in `celltype_new1` (= Celltype_Tier2 of Supplementary Data 4, 15 types).
- `tools/supp_tables_to_csv.py` converts the paper's supplementary Excel tables to CSV (used only to cross-check spot counts and MISTy results).
