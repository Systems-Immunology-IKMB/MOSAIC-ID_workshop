# Day 3 – Single-cell transcriptomics with a sex & gender lens

| File | Content |
|---|---|
| `D3_singleCell_Sex_2026.Rmd` | **Part 1**: intro (sex vs gender, COVID-19), loading 6 samples, QC, sex check (*XIST*/chrY), integration, clustering, annotation. Saves `output/D3_annotated_6samples.rds`. |
| `D3_singleCell_Sex_downstream_2026.Rmd` | **Part 2**: sex-sensitive downstream analyses: differential abundance (propeller), per-cell vs pseudobulk DE (DESeq2) with sex in the design, stratified and interaction contrasts, confounding demo, ORA + GSEA by sex, reporting checklist. |
| `sample_sheet.csv` | One row per sample. **Fill in the `FILL_ME` fields before running.** |

## Data layout

```
D3_Singlecell_Sex/
├── sample_sheet.csv
├── input/                       # not in git (patient data)
│   ├── <folder of sample 1>/raw_feature_bc_matrix/   (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
│   └── ...
└── output/                      # created by the notebooks
```

`sample_sheet.csv` columns:

- `sample_id`: short name used in all plots
- `folder`: the sub-folder in `input/` (e.g. `J09823_female_012_TA`), or the full path to the sample folder
- `donor`: the **person**. Two samples from the same person must have the same donor ID; Part 2 then adds the donor-blocking exercises automatically.

  Our design: 4 donors, one per sex × disease group. Each COVID-19 donor (one female, one male) was sampled at two time points (T1, T2). Healthy controls have one sample each. Donor IDs and time points are pre-filled; only the folder names are missing.
- `disease`: `healthy` / `COVID19`
- `sex`: `female` / `male`
- `timepoint`

References (Part 1 only, optional) are read from `params$ref_dir` (default `~/References`): `PBMCs_SeuratV4.rda`, `Housekeeping_Hs.csv`, `PBMCs_marker.csv`, `Reference_dmap.RDS`, `Reference_hpca.RDS`. If a file is missing, that step is skipped, and annotation falls back to canonical marker scores.

## Running

Student version (exercises without solutions):

```r
rmarkdown::render("D3_singleCell_Sex_2026.Rmd")
rmarkdown::render("D3_singleCell_Sex_downstream_2026.Rmd")
```

Instructor version (with solutions):

```r
rmarkdown::render("D3_singleCell_Sex_downstream_2026.Rmd",
                  params = list(solutions = TRUE), output_file = "D3_part2_instructor.html")
```

Other parameters: `input_dir`, `ref_dir`, `out_dir`. In Part 2, `celltype_of_interest` overrides the automatic choice and `min_cells` sets the pseudobulk minimum. Students who don't finish Part 1 can start Part 2 from the instructor's `D3_annotated_6samples.rds`.

## Packages

Seurat **v5** is required (`IntegrateLayers`, `JoinLayers`, layers).

```r
install.packages(c("Seurat", "dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork",
                   "pheatmap", "RColorBrewer", "rmarkdown", "BiocManager"))
BiocManager::install(c("speckle", "limma", "edgeR", "DESeq2", "clusterProfiler", "enrichplot",
                       "org.Hs.eg.db", "glmGamPoi",
                       "SingleR", "MAST"))   # the last two are optional
```

## What changed compared to the 2022 version

- **4 → 6 samples**, driven by the sample sheet instead of hard-coded objects. Metadata is added in a loop. In 2022 it was added to one sample only, and `SplitObject(..., "Patient")` failed for the others.
- `Gender` → `sex`, with an intro on the difference and a genetic **sex check** of every sample.
- Seurat v5 syntax: `layer =` instead of `slot =`, `IntegrateLayers(RPCA)` instead of `FindIntegrationAnchors`, `JoinLayers`.
- XIST/TSIX and chrY genes are removed from the variable genes so that clusters don't split by sex.
- Bug fixes:
  - `percent.hk`: the 2022 code passed a per-cell count vector as `features`. It is now the % of UMIs from housekeeping genes. **Please re-check the `> 10` threshold.**
  - `ScaleData` ran before `FindVariableFeatures`.
  - SingleR used `as.matrix()` on the full matrix (>10 GB with 6 samples).
  - The hpca reference path was inconsistent.
  - The cluster→cell-type mapping was hard-coded twice for a different dataset. It is now data-driven (majority vote of Azimuth l2 labels, or marker scores), with an exercise to annotate by hand.
- The 2022 per-cluster `FindAllMarkers` + GO section is **replaced** by DE-based ORA and GSEA per sex in Part 2. The 2022 Wilcoxon/bimod/MAST comparison is kept as Exercise 3.1.
- Sex colours are orange/green (colour-blind safe, no pink/blue).

## Testing

Both notebooks were rendered end to end (student and instructor versions) on **simulated** data: 6 samples, ~1,600 cells each, a planted sex signal, more CD14 monocytes in COVID-19 and in males, and a stronger interferon response in females. This was run with Seurat 5.0.1, DESeq2 1.42, clusterProfiler 4.10 and speckle 1.x. The reference-mapping and SingleR chunks could **not** be tested, because the reference files are not available here.
