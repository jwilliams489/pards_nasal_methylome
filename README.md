# Nasal Methylome Analysis — PARDS (EPIC v2.0)

Companion analysis code for Williams et al., *Methylomic Analysis of Nasal Brushings
Poorly Discriminates Subgroups in Pediatric Acute Respiratory Distress Syndrome*.

The pipeline is **deconvolution-first**: cell-type composition is estimated from the
methylation data and removed before any clustering, PCA, or differential-methylation
analysis, so that apparent "subgroup" structure is not driven by cellular composition.

## Requirements

- R 4.5.0
- Bioconductor + CRAN packages: `tidyverse`, `sesame`, `sesameData`, `BiocParallel`,
  `minfi`, `IlluminaHumanMethylationEPICv2manifest`,
  `IlluminaHumanMethylationEPICv2anno.20a1.hg38`, `limma`, `matrixStats`, `EpiDISH`,
  `cluster`, `clusterProfiler`, `msigdbr`, `AnnotationDbi`, `org.Hs.eg.db`,
  `GenomicRanges`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, `pheatmap`, `RColorBrewer`,
  `scatterplot3d`, `ggrepel`, `patchwork`, `wateRmelon`
- Optional (for formatted Word tables): `flextable`, `officer`

The script fetches the EPICv2 manifest and sesameData resources from the Zhou
InfiniumAnnotation mirror (`SESAMEDATA_USE_ALT = TRUE`), so an internet connection is
required on first run.

## Inputs

1. **IDAT files** — one `_Grn.idat` / `_Red.idat` pair per sample, in `Combined IDATs/`.
2. **Sample sheet** — `unified_samplesheet.csv` (one row per sample) with at least:
   `Sample_Name` (Sentrix `barcode_position`), `Subject_ID`, `Sentrix_ID`,
   `Sentrix_Position`, `Sample_Group` (ARDS/Control), `Sex`, `Age`, `Outcome`,
   plus the clinical fields used in Section 13 (e.g. PELOD, ventilator-free days,
   `Principal_Comorbidity`, `Infectious_Agent`).

## Running

1. Open `Methyl_Revis_Master.R` and set `PROJECT_ROOT` (top of the script) to the
   directory that contains `Combined IDATs/` and the sample sheet.
2. Source the script:

   ```r
   source("Methyl_Revis_Master.R")
   ```

Output directories (`figures_v2/`, `tables_v2/`, `logs_v2/`, `analysis_output_v2/`)
are created automatically under `PROJECT_ROOT`. A timestamped log and a full
`sessionInfo()` are written for provenance.

## Outputs

- **Figures** — unadjusted and adjusted PCA, hierarchical clustering, promoter-DMC
  pathway enrichment, and the supplementary QC and CellDMC leave-one-out panels.
- **Tables** — per-sample QC metrics, cell-type composition, silhouette diagnostics,
  clinical characteristics, adjusted-group DMC counts, and the CellDMC validation table.
- **`PARDS_methylome_full_analysis.rds`** — a single bundle of the key analysis objects
  (matrices, PCA, clustering, DMC results, CellDMC summaries) plus `sessionInfo()`.

## Notes

- The CellDMC permutation null (Section 15c) reruns CellDMC genome-wide `B_PERM` times
  and is the slow step. It defaults to `B_PERM = 100`; raise it (and `PERM_CORES`) for a
  reportable null.
- Raw IDATs and the processed beta matrix are deposited in NCBI GEO
  (accession: *to be added on acceptance*). Clinical metadata contain protected health
  information and are not shared publicly.
