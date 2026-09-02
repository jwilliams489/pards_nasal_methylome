# Nasal methylome analysis in pediatric ARDS (EPIC v2.0)

Companion analysis code for:

> Williams JG, Elmadany N, Joshi R, Jones R, Gregor N, Lahni P, Standage SW, Varisco BM.
> **Methylomic Analysis of Nasal Brushings in Pediatric Acute Respiratory Distress Syndrome: A Pilot Study.**

| | |
|---|---|
| Data | NCBI GEO [**GSE337899**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337899) |
| Archived code | [10.5281/zenodo.21245950](https://doi.org/10.5281/zenodo.21245950) |
| Platform | Illumina Infinium MethylationEPIC v2.0 |
| Specimens | 24 nasal brushings (21 PARDS, 3 control) |
| Analysis | R 4.5.0 |
| License | MIT (see `LICENSE`) |

The pipeline is deconvolution-first: cell-type composition is estimated from the
methylation data and removed before any PCA, clustering, or differential-methylation
analysis. The central finding is that composition accounts for more than 90% of the
variance in the unadjusted nasal methylome, and that no robust subgroup structure
survives adjustment for it at this sample size.

---

## Repository contents

| Path | What it is |
|---|---|
| `Methyl_PARDS_pipeline.R` | The complete analysis, raw iDATs to every reported figure and table. Runs start to finish. |
| `make_supplemental_figure_1.py` | Generates Supplemental Figure 1 (enrollment flow chart), which is drawn rather than computed. |
| `templates/unified_samplesheet_template.csv` | Sample sheet template: the exact header, plus two synthetic example rows to be replaced. |
| `LICENSE` | MIT. |

`Methyl_PARDS_pipeline.R` supersedes `Methyl_Revis_Master_Pub.R` from earlier commits.
It merges the original processing pipeline and the revision analyses into one file
and orders its sections to follow the Results.

---

## Getting the data

1. Download the raw iDATs from GEO series **GSE337899**. There are 48 files: two
   channels (`_Grn.idat`, `_Red.idat`) for each of 24 specimens. Both plain and
   gzipped (`.idat.gz`) files are read.
2. Put them under `<PROJECT_ROOT>/Combined IDATs/`. They may sit flat in that
   directory or stay nested in per-slide subdirectories — the script searches
   recursively.
3. Put `unified_samplesheet.csv` in the same directory (see below).

**File naming.** Each iDAT basename must correspond to a row in the sample sheet
by one of two keys, checked in this order:

| Key | Example basename |
|---|---|
| `Sample_Name` | `ARDS_002_Grn.idat` |
| `Sentrix_ID` + `_` + `Sentrix_Position` | `206891110001_R01C01_Grn.idat` |

A leading GEO sample-accession prefix is stripped before matching, so files
downloaded from GSE337899 unchanged — `GSM12345678_206891110001_R01C01_Grn.idat.gz`
— resolve on the Sentrix key without renaming. The script reports which key it
matched on, and if none resolves every specimen it names the specimens at fault
and shows example basenames it found on disk.

Expected layout:

```
PROJECT_ROOT/
└── Combined IDATs/
    ├── unified_samplesheet.csv
    ├── GSM12345678_206891110001_R01C01_Grn.idat.gz
    ├── GSM12345678_206891110001_R01C01_Red.idat.gz
    └── ...
```

Everything else — `figures_v2/`, `tables_v2/`, `logs_v2/`, `analysis_output_v2/`
— is created by the script after it has validated its inputs.

---

## Sample sheet

One row per specimen. `templates/unified_samplesheet_template.csv` carries the
exact header; replace its two `EXAMPLE_` rows with your own.

**What is and is not public.** The array data and the technical metadata in
GSE337899 — specimen identifiers, group, Sentrix barcode and position — are
enough to reproduce quality control, probe filtering, deconvolution, composition
adjustment, PCA, and clustering: Figures 1–5, Tables 2, 4, 5 and 7, and
Supplemental Tables 2–8. The clinical variables contain protected health
information and are not deposited. The analyses that depend on them — Table 1,
Table 3, and the immunocompromised contrast in Table 6 and Figure 6 — need those
fields, which are available from the corresponding author under an appropriate
data use agreement.

**Required — the script stops in seconds if any is absent:**

| Column | Contents |
|---|---|
| `Sample_Name` | Specimen key, e.g. `ARDS_002`. One of the two iDAT matching keys. |
| `Subject_ID` | Subject key. Equal to `Sample_Name` in this cohort. |
| `Sample_Group` | `ARDS` or `Control`. |
| `Sex` | `M`, `F`, or blank. Must be read as text, so do not leave the whole column empty. |
| `Age` | Years. Read by exact name for the covariate-adjusted CellDMC model in PART V. |
| `Principal_Comorbidity` | Free text. Read by exact name; the immunocompromised contrast in PART V is derived from it by keyword match. |

**Resolved by alias.** Each of these is looked up through a list of accepted
spellings, so any one will do. A group with no match drops that row from Table 1
or that test from Table 3, and the script warns rather than stopping.

| Variable | Accepted column names |
|---|---|
| Race | `Race`, `race` |
| Acute condition | `Acute_Condition`, `Acute_Dx`, `Acute_Diagnosis`, `Admission_Diagnosis` |
| PARDS severity | `PARDS_Severity`, `PARDS_Category`, `Highest_PARDS_Category`, `PARDS_severity`, `Severity`, `OI_Category` |
| PELOD-2 | `PELOD2`, `PELOD_2`, `PELOD`, `Highest_PELOD2`, `PELOD2_max` |
| Ventilator-free days | `VFD`, `Ventilator_Free_Days`, `VFD_28`, `Vent_Free_Days` |
| Outcome | `Outcome`, `Mortality`, `Survival`, `Vital_Status` |
| Array slide | `Sentrix_ID`, `Slide`, `SentrixID`, `Sentrix_Barcode` |
| Array position | `Sentrix_Position`, `Position`, `Sentrix_Pos` |

`Sentrix_ID` and `Sentrix_Position` are both a technical covariate in Table 3 and
the fallback iDAT matching key, so in practice include them under those names.

Specimens shared with our previous study (Williams et al., *Respir Res*
2022;23:181) are hard-coded as `DUPLICATE_SUBJECTS` in section I.1b rather than
read from the sheet, so the new-PARDS-only cohort is reproducible without it.

---

## Requirements

**R 4.5.0.** Every run writes the full package manifest to
`logs_v2/sessionInfo.txt`; the versions behind the published results are the ones
in that file from the archived run.

```r
install.packages(c("tidyverse", "matrixStats", "cluster", "ggrepel",
                   "patchwork", "vegan", "dendextend"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("sesame", "sesameData", "BiocParallel", "minfi", "limma",
                       "EpiDISH", "wateRmelon",
                       "IlluminaHumanMethylationEPICv2manifest",
                       "IlluminaHumanMethylationEPICv2anno.20a1.hg38"))

# Optional. Without these the script writes tables as CSV only; with them it
# also writes formatted .docx tables.
install.packages(c("flextable", "officer"))
```

The script checks every package before any analysis runs and reports all missing
ones at once with the install call, rather than failing at the first `library()`.

Earlier versions of this repository also required `clusterProfiler`, `msigdbr`,
`AnnotationDbi`, `org.Hs.eg.db`, `GenomicRanges`,
`TxDb.Hsapiens.UCSC.hg38.knownGene`, `pheatmap`, `scatterplot3d` and
`RColorBrewer`. Those supported a pathway-enrichment analysis of the clusters
that was removed from the manuscript during revision. None is called anywhere in
the current script and none is needed.

**Network access is required on first run.** `sesameAnno_download()` fetches
`EPICv2.hg38.manifest.tsv.gz` from the Zhou InfiniumAnnotation host to build the
iDAT address file, which pins the manifest version used for probe decoding. The
script also sets `options(SESAMEDATA_USE_ALT = TRUE)`, which routes sesameData
lookups to that same mirror; as of mid-2026 the EPIC v2 objects on ExperimentHub
redirect to an archived S3 prefix and cannot be retrieved. That option must stay
`TRUE` for the whole session.

---

## Running it

`PROJECT_ROOT` is the only path that needs to change. Set it in one of three ways:

```bash
# 1. environment variable, no edit to the file
PARDS_PROJECT_ROOT=/path/to/data Rscript Methyl_PARDS_pipeline.R
```

```r
# 2. from an R session
Sys.setenv(PARDS_PROJECT_ROOT = "/path/to/data")
source("Methyl_PARDS_pipeline.R")
```

```r
# 3. edit the CONFIGURATION block near the top of the script
PROJECT_ROOT <- "/path/to/data"
```

It defaults to `"."`, so running from inside a correctly laid-out directory works
with no configuration at all. Inputs are validated before any output directory is
created, so a run started in the wrong place fails cleanly without leaving
folders behind.

**Runtime.** `openSesame` on 24 EPIC v2 arrays takes roughly 20 minutes.
Everything except the permutation null runs in well under an hour.

**The permutation null (PART V) is the slow step.** One genome-wide CellDMC fit
takes about 8 minutes, so the full `B_PERM <- 1000` is roughly 5.5 days
single-threaded and is not a laptop job. Three ways through it:

- If `analysis_output_v2/celldmc_IC_permutation_null_B1000_full.rds` is present,
  it is loaded automatically and the permutation is skipped. This is checked
  before `RUN_PERMUTATION` is consulted, so a completed cluster run is picked up
  with no edit.
- Run it on a cluster. The loop checkpoints to disk every `CHUNK` (50)
  iterations and resumes where it stopped, so it survives job time limits.
- Set `RUN_PERMUTATION <- FALSE` to skip it. Table 6 is then produced without the
  empirical p-values.

---

## Outputs, and where each appears in the manuscript

Figures are written as both PDF (vector, submitted) and 600 dpi PNG.

| Manuscript item | Script section | Output file |
|---|---|---|
| Figure 1 | II | `figures_v2/Figure1_PCA_unadjusted.pdf` |
| Figure 2 | IV | `figures_v2/Figure2_PCA_adjusted.pdf` |
| Figure 3 | IV | `figures_v2/Figure3_dendrogram_adjusted.pdf` |
| Figure 4 | IV | `figures_v2/Figure4_PARDS_only_panels.pdf` |
| Figure 5 | IV | `figures_v2/Figure5_new_specimens_only_panels.pdf` |
| Figure 6 | V | `figures_v2/Figure6_CellDMC_leave_one_out.pdf` |
| Supplemental Figure 2 | I | `figures_v2/SuppFig2_QC_panel.pdf` |
| Table 1 | I | `tables_v2/Table1_cohort_descriptive.*` |
| Table 2 | III | `tables_v2/Table2_fraction_distributions.*` |
| Table 3 | III | `tables_v2/Table3_compositional_axis_tests.*` |
| Table 4 | III | `tables_v2/Table4_pc1_vs_composition.*` |
| Table 5 | IV | `tables_v2/Table5_silhouette_by_k_all_cohorts.*` |
| Table 6 | V | `tables_v2/Table6_CellDMC_immunocompromised.*` |
| Table 7 | VI | `tables_v2/Table7_adjustment_sensitivity.*` |
| Supplemental Table 2 | I | `tables_v2/SuppTable2_probe_filtering.*` |
| Supplemental Table 3 | IV | `tables_v2/SuppTable3_pca_variance_pre_post.*` |
| Supplemental Table 4 | III | `tables_v2/SuppTable4_reference_cpg_overlap.*` |
| Supplemental Tables 5–7 | III | `tables_v2/` per-specimen fractions, one per cohort |
| Supplemental Table 8 | III | `tables_v2/SuppTable8_deconvolution_goodness_of_fit.*` |

**Supplemental Figure 1** (enrollment flow chart) is drawn, not computed. Generate
it with `python3 make_supplemental_figure_1.py [outdir]` — requires `matplotlib`.
The enrollment counts are collected at the top of that file so they can be checked
against the Results without reading the layout code.

**Supplemental Table 1** (specimen log with duplicate status) is maintained by hand
and is not produced here. Its duplicate column is the source of
`DUPLICATE_SUBJECTS`.

`tables_v2/sensitivity/` holds supporting output reported in the response to
reviewers rather than in the manuscript: the PARDS-only covariate disposition,
per-specimen silhouette widths, per-probe variance explained, the stage-one
adjustment comparison, and the PERMANOVA against the 2022 transcriptomic
subgroups. PART VII is that material plus `sessionInfo()`.

---

## Analysis summary

| Part | What it does | Results section |
|---|---|---|
| I | QC, probe filtering, replicate collapse, deconvolution, composition adjustment, definition of the three analytic cohorts | Specimen Selection and Quality Control |
| II | Unadjusted PCA, all three cohorts | Dimensional Reduction in All Cohorts |
| III | Reference overlap, fraction distributions, model fit, compositional axis, PC1 correlations | Cell Type Deconvolution in All Cohorts |
| IV | Composition adjustment, PCA, hierarchical clustering | Adjustment of the Methylome for Cell Type |
| V | CellDMC, with leave-one-out and permutation validation | Cell-Type Specific Differentially Methylated Cytosines |
| VI | Six adjustment strategies compared in each cohort | Sensitivity Analysis of Different Cell Type Deconvolution Strategies |
| VII | Response-to-reviewer analyses and `sessionInfo()` | — |

Three nested analytic sets are carried throughout: the **full cohort** (24
specimens), **PARDS only** (21), and **new PARDS only** (13, excluding specimens
reused from the 2022 study). Deconvolution, composition adjustment, and the top
5% variance filter are re-derived independently within each. Only the probe set
is shared, because probe filtering is performed once across all specimens.

**Key parameters:**

| Constant | Value | Set in | Meaning |
|---|---|---|---|
| `VAR_THRESHOLD` | `1e-4` | I.1b | Cell fractions with variance below this are dropped from the covariate set |
| `B_PERM` | `1000` | I.1b | CellDMC permutation null |
| `N_PERM` | `9999` | I.1b | PERMANOVA permutations (PART VII) |
| `DETP_SAMPLE_MAX` | `0.05` | I.1 | Max mean pOOBAH detection p per specimen |
| `BSCON_MIN` | `80` | I.1 | Min bisulfite conversion, percent |
| `LOG2INT_MIN` | `10.5` | I.1 | Min mean log2 signal intensity |
| `IDENT_MATCH_MIN` | `0.90` | I.1 | Min rs-probe genotype concordance for the specimen-identity check |
| seed | `20260101` | I.1 | Set once at the top; the PERMANOVA resets it, and the CellDMC null seeds from its checkpoint so a resumed run reproduces |

Quality control flags a specimen on any of detection p, bisulfite conversion,
signal intensity, or the genotype-identity check. All 24 specimens passed every
check. The sesame GCT score is also computed but never contributes: it is not
returned for the EPIC v2 platform, so no specimen was ever evaluated against
`GCT_MAX`.

---

## Citation

Please cite the manuscript and the GEO accession GSE337899. For the code
specifically, cite the archived release: [10.5281/zenodo.21245950](https://doi.org/10.5281/zenodo.21245950).
