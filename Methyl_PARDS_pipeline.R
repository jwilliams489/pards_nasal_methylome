# =============================================================================
# Nasal methylome analysis in pediatric ARDS: full analysis pipeline
# Infinium MethylationEPIC v2.0
#
# Companion script for: Williams et al., "Methylomic Analysis of Nasal
# Brushings in Pediatric Acute Respiratory Distress Syndrome: A Pilot Study".
#
# This script merges the original processing pipeline (Methyl_Revis_Master.R)
# with the revision analyses (Methyl_Revis_Master2.R) into a single file that
# runs start to finish, from raw iDATs to every table and figure reported in
# the manuscript. Sections are ordered to follow the Results.
#
#   PART I    Processing: QC, filtering, replicate collapse, deconvolution,
#             composition adjustment, and definition of the three analytic
#             cohorts (Results: Specimen Selection and Quality Control)
#   PART II   Unadjusted PCA in all three cohorts
#             (Results: Dimensional Reduction in All Cohorts)
#   PART III  Cell-type deconvolution: reference overlap, fraction
#             distributions, model fit, compositional axis, PC1 correlations
#             (Results: Cell Type Deconvolution in All Cohorts)
#   PART IV   Composition adjustment, PCA, and clustering in all three cohorts
#             (Results: Adjustment of the Methylome for Cell Type)
#   PART V    Cell-type-specific differential methylation (CellDMC)
#             (Results: Cell-Type Specific Differentially Methylated Cytosines)
#   PART VI   Sensitivity to the deconvolution and adjustment strategy
#             (Results: Sensitivity Analysis of Different Cell Type
#             Deconvolution Strategies)
#   PART VII  Analyses reported in the response to reviewers only, and
#             session information
#
# Manuscript outputs produced, by section:
#   Figure 1  PART II    Unadjusted methylome PCA, full cohort
#   Figure 2  PART IV    PCA of the adjusted methylome, full cohort
#   Figure 3  PART IV    Hierarchical clustering of the adjusted methylome
#   Figure 4  PART IV    PCA and clustering, PARDS only
#   Figure 5  PART IV    PCA and clustering, new PARDS only
#   Figure 6  PART V     CellDMC leave-one-out
#   Table 1   PART I     Cohort characteristics
#   Table 2   PART III   Fraction distributions and boundary values
#   Table 3   PART III   Compositional log ratio vs clinical and technical
#   Table 4   PART III   Unadjusted PC1 vs cell fractions
#   Table 5   PART IV    Silhouette width and cluster sizes, all cohorts
#   Table 6   PART V     CellDMC counts and validation analyses
#   Table 7   PART VI    Adjustment strategy sensitivity
#   Supplemental Table 2  PART I    Probe filtering accounting
#   Supplemental Table 3  PART IV   PCA variance before and after adjustment
#   Supplemental Table 4  PART III  Reference CpGs retained
#   Supplemental Tables 5 to 7  PART III  Per-specimen fractions, each cohort
#   Supplemental Table 8  PART III  Deconvolution goodness of fit
#   Supplemental Figure 2 PART I    Per-sample array quality control
#   (Supplemental Figure 1, the enrollment flow chart, and Supplemental Table 1,
#    the specimen log, are maintained by hand and are not produced here. The
#    duplicate-status column of that log is the source of DUPLICATE_SUBJECTS in
#    section I.1b, which defines the new-PARDS-only cohort.)
#
# ORDERING NOTE. The Results report cell-type deconvolution after the
# unadjusted PCA, but the final probe set depends on the composition
# adjustment: probes whose residuals are non-finite after removeBatchEffect are
# dropped from all matrices so that every analysis runs on the same 298,445
# probes. Deconvolution and the full-cohort adjustment are therefore ESTIMATED
# at the end of PART I, and everything about them is REPORTED later, in
# manuscript order. No other section departs from the order of the Results.
#
# RUNTIME. openSesame on 24 EPIC v2 arrays takes roughly 20 minutes. The
# CellDMC permutation null (PART V) is the slow step: one genome-wide fit takes
# about 8 minutes, so B = 1000 permutations is not feasible on a laptop. Set
# RUN_PERMUTATION = FALSE to skip it, or run it on a cluster with the scripts
# in hpc/ (about 40 minutes across 20 array tasks). The local loop checkpoints
# every CHUNK iterations and resumes. A completed cluster run is picked up
# automatically, before RUN_PERMUTATION is consulted, if
# analysis_output_v2/celldmc_IC_permutation_null_B1000_full.rds is present.
# See PART V for details.
#
# DATA. Raw iDAT files and the per-specimen metadata are deposited in NCBI GEO
# under accession GSE337899. Code is archived at
# https://github.com/jwilliams489/pards_nasal_methylome. See README.md for the
# sample-sheet contract, package versions, and the output-to-manuscript map.
#
# REPRODUCIBILITY. This script is self-contained: every path is set in the
# CONFIGURATION block below, every required package is checked before any
# analysis runs, and every output directory is created if absent. The only
# edit needed for a fresh environment is PROJECT_ROOT (or the PARDS_PROJECT_ROOT
# environment variable).
#
# Author: James G. Williams, MD
# =============================================================================



# =============================================================================
# PART I. PROCESSING AND COHORT DEFINITION
# =============================================================================


# -----------------------------------------------------------------------------
# I.1 Setup, logging, and packages
# -----------------------------------------------------------------------------

# Every package this script actually calls, grouped by source. Earlier versions
# of this file also attached clusterProfiler, msigdbr, AnnotationDbi,
# org.Hs.eg.db, GenomicRanges, TxDb.Hsapiens.UCSC.hg38.knownGene, pheatmap,
# scatterplot3d and RColorBrewer. Those supported the pathway-enrichment
# analysis of the clusters, which was removed from the manuscript during
# revision; none of them is called anywhere below, and several are large, slow
# installs. They are deliberately no longer required.
CRAN_PKGS <- c("tidyverse", "matrixStats", "cluster", "ggrepel", "patchwork",
               "vegan", "dendextend")
BIOC_PKGS <- c("sesame", "sesameData", "BiocParallel", "minfi", "limma",
               "EpiDISH", "wateRmelon",
               "IlluminaHumanMethylationEPICv2manifest",
               "IlluminaHumanMethylationEPICv2anno.20a1.hg38")
# Optional: formatted .docx tables. The script writes CSV either way.
OPTIONAL_PKGS <- c("flextable", "officer")

# Report EVERY missing package at once, with the exact install call, rather
# than dying at the first library() and making the reader iterate.
.missing <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
.miss_cran <- .missing(CRAN_PKGS)
.miss_bioc <- .missing(BIOC_PKGS)
if (length(.miss_cran) || length(.miss_bioc)) {
  msg <- "Required packages are not installed.\n"
  if (length(.miss_cran))
    msg <- paste0(msg, "\n  install.packages(c(",
                  paste(sprintf('"%s"', .miss_cran), collapse = ", "), "))\n")
  if (length(.miss_bioc))
    msg <- paste0(msg,
      "\n  if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
      "install.packages(\"BiocManager\")\n",
      "  BiocManager::install(c(",
      paste(sprintf('"%s"', .miss_bioc), collapse = ", "), "))\n")
  stop(msg, call. = FALSE)
}
.miss_opt <- .missing(OPTIONAL_PKGS)
if (length(.miss_opt))
  message("Optional packages absent (tables will be written as CSV only): ",
          paste(.miss_opt, collapse = ", "),
          "\n  install.packages(c(",
          paste(sprintf('"%s"', .miss_opt), collapse = ", "), "))")

suppressPackageStartupMessages({
  library(tidyverse)
  library(sesame)
  library(sesameData)
  library(BiocParallel)  # SerialParam(): required by openSesame BPPARAM
  library(minfi)
  library(IlluminaHumanMethylationEPICv2manifest)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  library(limma)
  library(matrixStats)
  library(EpiDISH)
  library(cluster)
  library(ggrepel)       # point labels in QC figures
  library(patchwork)     # multi-panel supplementary QC figure
  library(wateRmelon)    # bscon(): bisulfite conversion % from control probes
})

# -----------------------------------------------------------------------------
# Namespace conflicts: pin the tidyverse verbs this script uses
# -----------------------------------------------------------------------------
# tidyverse is attached FIRST, so every Bioconductor package loaded after it
# that re-exports a tidyverse verb name wins the binding. The one that bites
# here is matrixStats::count(), which masks dplyr::count() and fails on a data
# frame with "Argument 'x' is not a vector: list". S4Vectors::rename() and
# IRanges::slice() are the same hazard a few hundred lines later. Rebinding is
# explicit and adds no dependency; conflicted::conflict_prefer_all("dplyr")
# would do the same job if you would rather declare it.
count  <- dplyr::count
rename <- dplyr::rename
slice  <- dplyr::slice
filter <- dplyr::filter
desc   <- dplyr::desc
first  <- dplyr::first
reduce <- purrr::reduce

options(ExperimentHub.ASK = FALSE)

# -----------------------------------------------------------------------------
# Route sesameData lookups to the alternate host (NOT the frozen ExperimentHub)
# -----------------------------------------------------------------------------
# As of mid-2026 the EPICv2 sesameData objects on ExperimentHub redirect to an
# AWS S3 prefix that has been moved to Glacier DEEP_ARCHIVE, so any internal
# sesameDataGet() of them throws "either not found or needs to be cached".
# Critically, several lookups are buried inside sesame functions we cannot reach
# with arguments -- most importantly pOOBAH(), whose nonuniqMask() pulls
# KYCG.EPICv2.Mask.20230314 to exclude multi-mapping probes from its detection
# null (confirmed via traceback). Setting SESAMEDATA_USE_ALT = TRUE makes
# .sesameDataGet resolve titles from sesame's alternate (Zhou-hosted) mirror and
# cache them in-session, bypassing the archived bucket. NOTE: this is a global,
# session-wide reroute and must stay TRUE for the whole run -- toggling it off
# makes .sesameDataGet look up cacheEnv by EH-id rather than title, which would
# orphan any object the fallback stored under its title string.
options(SESAMEDATA_USE_ALT = TRUE)

set.seed(20260101)

# -----------------------------------------------------------------------------
# Configuration (edit PROJECT_ROOT to your environment)
# -----------------------------------------------------------------------------
# PROJECT_ROOT is the ONLY path that must be changed to run this script
# elsewhere. Everything below is derived from it. It may also be supplied
# without editing the file at all:
#
#   Sys.setenv(PARDS_PROJECT_ROOT = "/path/to/data"); source("Methyl_PARDS_pipeline.R")
#   PARDS_PROJECT_ROOT=/path/to/data Rscript Methyl_PARDS_pipeline.R
#
# The default, ".", runs the script from inside a directory laid out as
# described in README.md, which is what a fresh clone plus a GEO download gives.
PROJECT_ROOT <- Sys.getenv("PARDS_PROJECT_ROOT", unset = ".")
# The 48 iDAT files (24 specimens, two channels) sit flat in IDAT_DIR. They
# were consolidated there from the per-slide OneDrive delivery folders; section
# I.3 searches IDAT_DIR recursively and matches on file basename, so either
# layout works and no path edit is needed if the arrays are re-nested. The
# sample sheet path is set separately so the two can be moved independently.
IDAT_DIR     <- file.path(PROJECT_ROOT, "Combined IDATs")
SHEET_DIR    <- file.path(PROJECT_ROOT, "Combined IDATs")
# The unified sample sheet: one row per specimen, carrying the Sentrix chip
# columns and the clinical fields. The exact contract, including the accepted
# alternative spellings for each clinical variable, is documented in README.md
# and in the SHEET_* checks immediately below; a header-only template ships in
# templates/unified_samplesheet_template.csv. The same information is deposited
# as the per-sample characteristics of GEO series GSE337899, so the sheet can
# be rebuilt from the public record.
SAMPLE_SHEET <- file.path(SHEET_DIR, "unified_samplesheet.csv")
CLINICAL     <- NULL  # clinical data already merged into the unified sheet
FIG_DIR      <- file.path(PROJECT_ROOT, "figures_v2")
TBL_DIR      <- file.path(PROJECT_ROOT, "tables_v2")
LOG_DIR      <- file.path(PROJECT_ROOT, "logs_v2")
RDS_DIR      <- file.path(PROJECT_ROOT, "analysis_output_v2")
# Validate the INPUTS before creating any output directory. Creating them first
# means a run started in the wrong directory -- easy now that PROJECT_ROOT
# defaults to "." -- scatters four empty folders into the working tree and only
# then reports the real problem.
for (p_ in c(IDAT_DIR, SHEET_DIR))
  if (!dir.exists(path.expand(p_)))
    stop("Directory does not exist: ", p_,
         "\n  Set PROJECT_ROOT in the block above, or run with",
         "\n    PARDS_PROJECT_ROOT=/path/to/data Rscript Methyl_PARDS_pipeline.R",
         "\n  Expected layout is documented in README.md.")
if (!file.exists(path.expand(SAMPLE_SHEET)))
  stop("Sample sheet not found: ", SAMPLE_SHEET,
       "\n  See README.md for the required columns; a template with the exact",
       "\n  header is provided at templates/unified_samplesheet_template.csv.")

# The sample-sheet contract, checked here -- before 20 minutes of openSesame --
# so that a bad sheet fails in seconds rather than deep inside PART I.
#
# Two tiers, matching how the script actually reads the sheet. SHEET_REQUIRED
# columns are addressed by exact name somewhere below and their absence is
# fatal -- note that Age and Principal_Comorbidity belong here even though
# Table 1 resolves them through aliases, because PART V reads ss_s3$Age and
# sample_sheet$Principal_Comorbidity literally. The variables in SHEET_ALIASES
# are only ever reached through the alias lists that the Table 1 builder
# (VAR_MAP) and the technical-covariate picker (pick1) use, so any one spelling
# will do; a missing group drops that row from Table 1 or that test from
# Table 3, so it warns rather than stops.
SHEET_REQUIRED <- c("Sample_Name",           # matched to the iDAT basename
                    "Subject_ID",
                    "Sample_Group",          # CASE_LABEL / CONTROL_LABEL
                    "Sex",                   # PART I reads sample_sheet$Sex
                    "Age",                   # PART V covariate model reads $Age
                    "Principal_Comorbidity") # source of the PART V IC contrast
SHEET_ALIASES <- list(
  `Race`                  = c("Race", "race"),
  `Acute condition`       = c("Acute_Condition", "Acute_Dx", "Acute_Diagnosis",
                              "Admission_Diagnosis"),
  `PARDS severity`        = c("PARDS_Severity", "PARDS_Category",
                              "Highest_PARDS_Category", "PARDS_severity",
                              "Severity", "OI_Category"),
  `PELOD-2`               = c("PELOD2", "PELOD_2", "PELOD", "Highest_PELOD2",
                              "PELOD2_max"),
  `Ventilator-free days`  = c("VFD", "Ventilator_Free_Days", "VFD_28",
                              "Vent_Free_Days"),
  `Outcome`               = c("Outcome", "Mortality", "Survival",
                              "Vital_Status"),
  `Array slide`           = c("Sentrix_ID", "Slide", "SentrixID",
                              "Sentrix_Barcode"),
  `Array position`        = c("Sentrix_Position", "Position", "Sentrix_Pos")
)
.hdr <- names(readr::read_csv(SAMPLE_SHEET, n_max = 0, show_col_types = FALSE))
.missing_req <- setdiff(SHEET_REQUIRED, .hdr)
if (length(.missing_req))
  stop("Sample sheet is missing required column(s): ",
       paste(.missing_req, collapse = ", "),
       "\n  File: ", SAMPLE_SHEET,
       "\n  See README.md, section \"Sample sheet\".")
.unmatched <- names(SHEET_ALIASES)[
  !vapply(SHEET_ALIASES, function(a) any(a %in% .hdr), logical(1))]
if (length(.unmatched))
  warning("Sample sheet has no column for: ", paste(.unmatched, collapse = ", "),
          ".\n  These rows will be dropped from Table 1 or Table 3. ",
          "Accepted spellings are listed in README.md.", call. = FALSE)
rm(.hdr, .missing_req, .unmatched, .missing, .miss_cran, .miss_bioc, .miss_opt)

# Inputs check out; now create the output tree. Unlike the loop this replaces,
# a failure here is reported at once rather than surfacing at the first ggsave.
for (d in c(FIG_DIR, TBL_DIR, LOG_DIR, RDS_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
if (!all(dir.exists(c(FIG_DIR, TBL_DIR, LOG_DIR, RDS_DIR))))
  stop("Could not create output directories under: ", PROJECT_ROOT,
       "\n  Check that the path exists and is writable.")

# -----------------------------------------------------------------------------
# Publication-ready table writer
# -----------------------------------------------------------------------------
# Always writes a CSV. If flextable + officer are installed, also writes a
# formatted Word (.docx) table with a numbered caption and footnote, ready to
# drop into the manuscript. Install once with:
#   install.packages(c("flextable", "officer"))
save_pub_table <- function(df, file_noext, caption = NULL, note = NULL) {
  readr::write_csv(df, paste0(file_noext, ".csv"))
  if (requireNamespace("flextable", quietly = TRUE) &&
      requireNamespace("officer", quietly = TRUE)) {
    ft <- flextable::flextable(df)
    ft <- flextable::theme_booktabs(ft)
    ft <- flextable::align(ft, part = "all", align = "center")
    ft <- flextable::align(ft, j = 1, part = "all", align = "left")
    ft <- flextable::font(ft, fontname = "Arial", part = "all")
    ft <- flextable::fontsize(ft, size = 9, part = "all")
    ft <- flextable::bold(ft, part = "header")
    ft <- flextable::padding(ft, padding = 3, part = "all")
    ft <- flextable::autofit(ft)
    doc <- officer::read_docx()
    if (!is.null(caption))
      doc <- officer::body_add_par(doc, caption, style = "heading 2")
    doc <- flextable::body_add_flextable(doc, ft)
    if (!is.null(note))
      doc <- officer::body_add_par(doc, note, style = "Normal")
    print(doc, target = paste0(file_noext, ".docx"))
    cat("  Wrote publication table:", basename(paste0(file_noext, ".docx")), "\n")
  } else {
    cat("  (flextable/officer not installed -> wrote CSV only for",
        basename(paste0(file_noext, ".csv")),
        "; install them for the formatted .docx)\n")
  }
  invisible(paste0(file_noext, ".csv"))
}

# QC thresholds for Tier-1 reviewer metrics (tune to your cohort before reporting)
DETP_SAMPLE_MAX <- 0.05    # max mean pOOBAH detection p per sample (operative exclusion)
GCT_MAX         <- 1.2     # bisulfite GCT: ~1.0 = complete, elevated = incomplete
LOG2INT_MIN     <- 10.5    # min mean log2 {M,U} intensity (minfi badSampleCutoff default)
BSCON_MIN       <- 80      # min bisulfite conversion % (wateRmelon::bscon)

# Genotype identity (EPIC rs/SNP probes) — used for swap / replicate checks
IDENT_MATCH_MIN <- 0.90    # min genotype concordance to call two samples the same person
RS_HET_LOW      <- 0.25    # rs beta < this -> homozygous (call 0)
RS_HET_HIGH     <- 0.75    # rs beta > this -> homozygous (call 2); between -> het (call 1)
RS_MIN_PROBES   <- 10      # minimum complete rs probes required to attempt the check

log_file <- file.path(LOG_DIR, paste0("methylation_v2_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)
cat("Methylation analysis (v2) started:", format(Sys.time()), "\n")
cat("R version:", R.version.string, "\n\n")

# -----------------------------------------------------------------------------
# EPICv2 address file: local download (belt-and-suspenders for iDAT decoding)
# -----------------------------------------------------------------------------
# The cold-archive problem (see SESAMEDATA_USE_ALT note above) is fully handled
# by routing sesameDataGet to the alternate host. We ALSO build the iDAT address
# file from the Zhou manifest TSV and pass it explicitly to openSesame() below,
# for two reasons: (1) it pins the manifest version used for decoding (good
# provenance: EPICv2.hg38.manifest.tsv.gz, Zhou InfiniumAnnotation), and (2) it
# was verified to decode these iDATs this session independent of the cache.
# Masking is NOT built here: pOOBAH's internal nonuniqMask and section 6 both
# resolve the real KYCG.EPICv2.Mask.20230314 via the alternate host.
cat("Downloading EPICv2 manifest from the Zhou InfiniumAnnotation host...\n")
epicv2_tsv  <- sesameAnno_download("EPICv2.hg38.manifest.tsv.gz")
epicv2_addr <- sesameAnno_buildAddressFile(epicv2_tsv)
cat("  Address file built from EPICv2.hg38.manifest.tsv.gz\n")



# -----------------------------------------------------------------------------
# I.1b Revision output directories, analysis constants, and shared helpers
# -----------------------------------------------------------------------------
# Sensitivity and supporting output is kept out of the primary figure and table
# directories so those hold only what the manuscript cites.
SENS_FIG_DIR <- file.path(FIG_DIR, "sensitivity")
SENS_TBL_DIR <- file.path(TBL_DIR, "sensitivity")
for (d in c(SENS_FIG_DIR, SENS_TBL_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Group labels as coded in the unified sample sheet.
CASE_LABEL    <- "ARDS"
CONTROL_LABEL <- "Control"

# Composition adjustment: any cell-type fraction whose variance across
# specimens falls below this threshold is dropped from the covariate set.
# Stated in the Methods; used identically in every cohort and every strategy.
VAR_THRESHOLD <- 1e-4

# Permutations for the PERMANOVA in PART VII (fast, one distance matrix) and
# for the CellDMC null in PART V (slow, one genome-wide fit per permutation).
N_PERM <- 9999
B_PERM <- 1000

# Specimens shared with Williams et al., Respir Res 2022;23:181, per
# Supplemental Table 1. Excluded in the new-PARDS-only cohort.
DUPLICATE_SUBJECTS <- c("ARDS_001", "ARDS_003", "ARDS_004", "ARDS_010",
                        "ARDS_015", "ARDS_018", "ARDS_019", "ARDS_021",
                        "Control_001", "Control_022")

# Shared plotting constants.
#
# DIAGNOSIS colours (Figures 1 and 2): orange = PARDS, blue = control.
FIG1_ORANGE <- "#D55E00"
FIG1_BLUE   <- "#56B4E9"

# NOTE ON CROSS-FIGURE READING: adj_assign, cl_pards and cl_new are three
# INDEPENDENT cutree() labellings of three different cohorts, so "Group 1" in
# Figure 3 is not the same set of subjects as "Group 1" in Figure 4 (the
# concordance table built in IV.5 exists precisely because they diverge).
# Sharing one palette across the three figures makes them look comparable, so
# the figure legends must state that group numbering is cohort-specific.
#
# CLUSTER colours (Figures 3, 4, 5). Deliberately NOT the diagnosis colours:
# FIG1_ORANGE and FIG1_BLUE mean "PARDS" and "control" in Figures 1 and 2, so
# reusing them for cluster membership invites the reader to read Group 1 as the
# PARDS group. Okabe-Ito, colourblind safe. To revert to the diagnosis palette
# set CLUSTER_PALETTE <- c(FIG1_ORANGE, FIG1_BLUE, "#009E73", "#E69F00").
CLUSTER_PALETTE <- c("#009E73", "#CC79A7", "#0072B2", "#E69F00",
                     "#000000", "#999999")

# Specimen labels in figures. Every table uses the Subject_ID form (ARDS_002),
# so figures use it too. Set TRUE for the compact ARDS002 form used previously.
STRIP_UNDERSCORE <- FALSE
fig_label <- function(x) if (STRIP_UNDERSCORE) gsub("_", "", x) else x

# Cluster legend labels, shared by Figures 3, 4 and 5 so their wording cannot
# drift apart. Takes a cutree() vector, returns a named label lookup.
cluster_labels <- function(cl) {
  sz <- table(cl)
  setNames(sprintf("Group %s (n = %d)", names(sz), as.integer(sz)), names(sz))
}

# Display names for the three-compartment CellDMC collapse, shared by Figure 6
# and Table 6. This covers ONLY that collapse; the stage-one cell-type naming
# used by the composition tables is a different grouping and is left alone.
COMPARTMENT_LABELS <- c(Epi = "Epithelial", Neutro = "Neutrophils",
                        Rest = "Others")

# Usable text width of the manuscript template (A4 portrait, 0.98 in margins).
# Set to 6.5 for US Letter with 1 in margins.
USABLE_WIDTH_IN <- 6.30

suppressPackageStartupMessages({
  library(vegan)       # adonis2(), betadisper(), vegdist() (PART VII)
  library(dendextend)  # as.ggdend() for the paneled dendrograms
})

# Adjusted Rand index, written out to avoid adding a dependency (mclust/fossil)
# for a single number. Compares two partitions of the same specimens after
# correcting for the agreement expected by chance.
adjusted_rand_index <- function(a, b) {
  tab <- table(a, b)
  n   <- sum(tab)
  sum_comb <- function(x) sum(choose(x, 2))
  idx      <- sum_comb(as.vector(tab))
  exp_idx  <- sum_comb(rowSums(tab)) * sum_comb(colSums(tab)) / choose(n, 2)
  max_idx  <- (sum_comb(rowSums(tab)) + sum_comb(colSums(tab))) / 2
  (idx - exp_idx) / (max_idx - exp_idx)
}

# Top 5% variance filter + PCA. One implementation, used for every PCA in the
# manuscript so the filter can never drift between sections.
pca_top5 <- function(mat) {
  pv    <- rowVars(mat, na.rm = TRUE)
  m_top <- mat[pv >= quantile(pv, 0.95, na.rm = TRUE), , drop = FALSE]
  pca   <- prcomp(t(m_top), center = TRUE, scale. = FALSE)
  list(m_top = m_top, pca = pca,
       var_exp = (pca$sdev^2) / sum(pca$sdev^2))
}

# Top 5% variance filter, PCA, Ward.D2 tree, silhouette over k. Used for every
# clustering analysis (full cohort, PARDS only, new PARDS only, and each
# adjustment strategy in PART VI).
cluster_diagnostics <- function(m, label, k_max = 6) {
  p       <- pca_top5(m)
  m_top   <- p$m_top
  var_exp <- p$var_exp
  cat("\n[", label, "] top 5% variable CpGs:", nrow(m_top), "\n")
  cat("[", label, "] variance explained (PC1 to PC5): ",
      paste0(round(var_exp[1:min(5, length(var_exp))] * 100, 1), "%",
             collapse = ", "), "\n", sep = "")

  d  <- dist(t(m_top), method = "euclidean")
  hc <- hclust(d, method = "ward.D2")

  k_top <- min(k_max, ncol(m) - 1)
  sil <- map_dfr(2:k_top, function(k) {
    cl <- cutree(hc, k = k)
    if (length(unique(cl)) < k) return(tibble(k = k, mean_silhouette = NA_real_))
    tibble(k = k, mean_silhouette = mean(silhouette(cl, d)[, 3]))
  })

  best_k   <- sil$k[which.max(sil$mean_silhouette)]
  best_sil <- max(sil$mean_silhouette, na.rm = TRUE)
  cat("[", label, "] mean silhouette by k:\n", sep = "")
  print(sil %>% mutate(mean_silhouette = round(mean_silhouette, 3)))
  cat("[", label, "] best k = ", best_k, " (width = ", round(best_sil, 3), ")\n",
      sep = "")
  cat("Interpretation: > 0.5 strong, 0.25-0.5 weak, < 0.25 essentially absent\n")

  list(m_top = m_top, pca = p$pca, var_exp = var_exp, dist = d, hc = hc,
       sil = sil, best_k = best_k, best_sil = best_sil, label = label)
}


# Supporting output. Same signature as save_pub_table(), but writes only the
# CSV: these analyses answer a reviewer comment or document an intermediate
# step, and are not manuscript tables. Keeping them as data rather than as
# formatted documents keeps every number auditable without filling the output
# directory with files that could be mistaken for submission material.
save_support_table <- function(df, file_noext, caption = NULL, note = NULL) {
  readr::write_csv(df, paste0(file_noext, ".csv"))
  cat("  Wrote supporting table:", basename(paste0(file_noext, ".csv")), "\n")
  invisible(paste0(file_noext, ".csv"))
}

# -----------------------------------------------------------------------------
# Grouped narrow-table writer (Tables 2 and 7)
# -----------------------------------------------------------------------------
# save_pub_table() calls flextable::autofit(), which sizes columns to their
# content and overruns the text block for tables that repeat a cohort label on
# every row. This writer instead promotes that label to a bold merged section
# row and sets explicit column widths under a fixed layout, so Word wraps
# inside the cells rather than widening the table. Widths default to
# proportional-to-content, rescaled to the usable text width; pass `widths` to
# override. The CSV keeps the ungrouped data frame, including the group column.
save_grouped_table <- function(df, group_col, file_noext, caption = NULL,
                               note = NULL, widths = NULL, font_size = 8,
                               page_width_in = USABLE_WIDTH_IN) {
  stopifnot(group_col %in% names(df))
  readr::write_csv(df, paste0(file_noext, ".csv"))

  data_cols <- setdiff(names(df), group_col)
  groups    <- unique(as.character(df[[group_col]]))

  body <- map_dfr(groups, function(g) {
    hdr <- as_tibble(setNames(as.list(c(g, rep("", length(data_cols) - 1L))),
                              data_cols)) %>% mutate(.row_type = "group")
    dat <- df %>%
      filter(as.character(.data[[group_col]]) == g) %>%
      dplyr::select(all_of(data_cols)) %>%
      mutate(across(everything(), as.character), .row_type = "data")
    bind_rows(hdr, dat)
  })
  row_type <- body$.row_type
  body     <- dplyr::select(body, -.row_type)

  if (!(requireNamespace("flextable", quietly = TRUE) &&
        requireNamespace("officer", quietly = TRUE))) {
    cat("  (flextable/officer not installed -> wrote CSV only for",
        basename(paste0(file_noext, ".csv")), ")\n")
    return(invisible(paste0(file_noext, ".csv")))
  }

  if (is.null(widths)) {
    w <- vapply(seq_along(data_cols), function(j)
      max(nchar(c(data_cols[j], body[[j]])), na.rm = TRUE), numeric(1))
    widths <- pmax(w, 6)
  }
  stopifnot(length(widths) == ncol(body))
  widths <- widths / sum(widths) * page_width_in

  grp_i <- which(row_type == "group")
  dat_i <- which(row_type == "data")
  rule  <- officer::fp_border(color = "black",  width = 0.75)
  hair  <- officer::fp_border(color = "grey60", width = 0.5)

  ft <- flextable::flextable(body)
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::font(ft, fontname = "Arial", part = "all")
  ft <- flextable::fontsize(ft, size = font_size, part = "all")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::align(ft, j = 1, align = "left", part = "all")
  ft <- flextable::align(ft, j = 2:ncol(body), align = "center", part = "all")
  # group rows: bold, merged across the full width, hairline above each block
  ft <- flextable::merge_h_range(ft, i = grp_i, j1 = 1, j2 = ncol(body),
                                 part = "body")
  ft <- flextable::bold(ft, i = grp_i, j = 1, part = "body")
  if (length(grp_i) > 1)
    ft <- flextable::hline(ft, i = grp_i[-1] - 1, border = hair, part = "body")
  ft <- flextable::padding(ft, i = dat_i, j = 1, padding.left = 12, part = "body")
  ft <- flextable::padding(ft, padding.top = 1.5, padding.bottom = 1.5,
                           part = "body")
  ft <- flextable::hline_top(ft, border = rule, part = "header")
  ft <- flextable::hline_bottom(ft, border = rule, part = "body")
  ft <- flextable::width(ft, j = seq_along(widths), width = widths)
  ft <- flextable::set_table_properties(ft, layout = "fixed")

  doc <- officer::read_docx()
  if (!is.null(caption))
    doc <- officer::body_add_par(doc, caption, style = "heading 2")
  doc <- flextable::body_add_flextable(doc, ft)
  if (!is.null(note))
    doc <- officer::body_add_par(doc, note, style = "Normal")
  print(doc, target = paste0(file_noext, ".docx"))
  cat("  Wrote publication table:", basename(paste0(file_noext, ".docx")),
      sprintf("(%.2f in wide)\n", sum(widths)))
  invisible(paste0(file_noext, ".csv"))
}


# -----------------------------------------------------------------------------
# I.2 Load sample sheet
# -----------------------------------------------------------------------------

# sesame's openSesame() handles EPIC v2.0 manifest properly. The Welsh 2023
# benchmark recommends "noob" + nonlinear dye-bias correction for EPIC v2.

sample_sheet <- read_csv(SAMPLE_SHEET, show_col_types = FALSE) %>%
  mutate(Sex = na_if(Sex, "NA"))

cat("Loaded unified sample sheet:", nrow(sample_sheet), "samples,",
    ncol(sample_sheet), "columns\n")
cat("Columns:", paste(names(sample_sheet), collapse = ", "), "\n")
cat("Sample_Group counts:\n")
print(table(sample_sheet$Sample_Group, useNA = "ifany"))



# -----------------------------------------------------------------------------
# I.3 Process iDATs with sesame openSesame (Noob + dye-bias)
# -----------------------------------------------------------------------------

# openSesame() returns a beta matrix with EPICv2 IDs (including the _TC21
# / _BC21 replicate suffixes). They are collapsed in section I.7.

# iDAT prefixes. The files are stored one directory per slide, and file names
# are matched case-insensitively, so this resolves whether the arrays sit flat
# in IDAT_DIR or nested. Every specimen is verified to have both channels
# before openSesame() runs.
GRN_PAT <- "_Grn\\.idat(\\.gz)?$"

grn_files <- list.files(IDAT_DIR, pattern = GRN_PAT, recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
if (length(grn_files) == 0)
  stop("No _Grn.idat files found under ", IDAT_DIR)

grn_base <- sub(GRN_PAT, "", basename(grn_files), ignore.case = TRUE)
if (any(duplicated(grn_base)))
  warning("Duplicate array basenames under ", IDAT_DIR, ": ",
          paste(unique(grn_base[duplicated(grn_base)]), collapse = ", "),
          ". The first match on disk is used for each.")

# Resolve each specimen to its array file. Our own iDATs are named for the
# specimen (ARDS_002_Grn.idat); files downloaded from GEO are named for the
# array and carry a sample-accession prefix
# (GSM12345678_206891110001_R01C01_Grn.idat.gz). Both are supported, so the
# script runs on an unmodified GEO download with no renaming: two sheet keys
# are tried against two spellings of each disk basename, and the first
# combination that resolves every specimen wins.
grn_keys <- list(
  `file name`      = grn_base,
  `GEO, GSM prefix stripped` = sub("^GSM[0-9]+_", "", grn_base)
)
sheet_keys <- list(`Sample_Name` = as.character(sample_sheet$Sample_Name))
if (all(c("Sentrix_ID", "Sentrix_Position") %in% names(sample_sheet)))
  sheet_keys[["Sentrix_ID_Sentrix_Position"]] <-
    paste0(sample_sheet$Sentrix_ID, "_", sample_sheet$Sentrix_Position)

hit <- NULL
for (sk in names(sheet_keys)) {
  for (gk in names(grn_keys)) {
    cand <- match(sheet_keys[[sk]], grn_keys[[gk]])
    if (!anyNA(cand)) {
      hit <- cand
      cat("Matched iDATs by ", sk, " against ", gk, ".\n", sep = "")
      break
    }
  }
  if (!is.null(hit)) break
}

if (is.null(hit)) {
  # Report against the sheet key that got closest, so the message names the
  # specimens actually at fault rather than all 24.
  best <- which.min(vapply(sheet_keys, function(s)
    min(vapply(grn_keys, function(g) sum(is.na(match(s, g))), numeric(1))),
    numeric(1)))
  bs   <- sheet_keys[[best]]
  bg   <- grn_keys[[which.min(vapply(grn_keys, function(g)
            sum(is.na(match(bs, g))), numeric(1)))]]
  missing <- bs[is.na(match(bs, bg))]
  stop("No Grn iDAT found for ", length(missing), " of ", nrow(sample_sheet),
       " specimens, matching on ", names(sheet_keys)[best], ": ",
       paste(head(missing, 5), collapse = ", "),
       if (length(missing) > 5) ", ..." else "",
       "\n  Searched recursively under: ", IDAT_DIR,
       "\n  Example basenames on disk: ",
       paste(head(grn_base, 3), collapse = ", "),
       "\n  Each iDAT basename must equal either the Sample_Name or the",
       "\n  Sentrix_ID_Sentrix_Position of its row in the sample sheet",
       "\n  (a leading GSM accession prefix is ignored). See README.md.")
}

idat_prefixes <- sub(GRN_PAT, "", grn_files[hit], ignore.case = TRUE)
names(idat_prefixes) <- sample_sheet$Sample_Name

no_red <- idat_prefixes[!file.exists(paste0(idat_prefixes, "_Red.idat")) &
                        !file.exists(paste0(idat_prefixes, "_Red.idat.gz"))]
if (length(no_red) > 0)
  stop("Grn present but Red missing for: ",
       paste(names(no_red), collapse = ", "))

cat("Located", length(idat_prefixes), "iDAT pairs under ", IDAT_DIR, "\n")
# Index-based, so this stays correct whichever key resolved the match above.
extra <- grn_base[-hit]
if (length(extra) > 0)
  cat("  Note:", length(extra), "arrays are present on disk but not in the",
      "sample sheet and are not processed.\n")

cat("Running sesame openSesame() on", length(idat_prefixes),
    "samples (Noob + dye-bias correction)...\n")

beta_raw <- openSesame(
  idat_prefixes,
  prep      = "QCDPB",     # Quality mask + Channel + Dye + Pval + Bg = Noob-like
  func      = getBetas,
  platform  = "EPICv2",    # explicit platform (decode via local address file)
  manifest  = epicv2_addr, # local address file (pins manifest version)
  BPPARAM   = SerialParam()
)

# openSesame returns a list of vectors when given a vector; coerce to matrix
if (is.list(beta_raw)) {
  common <- Reduce(intersect, lapply(beta_raw, names))
  beta_raw <- do.call(cbind, lapply(beta_raw, function(v) v[common]))
}
colnames(beta_raw) <- sample_sheet$Sample_Name

cat("Beta matrix:", dim(beta_raw), "\n")
saveRDS(beta_raw, file.path(RDS_DIR, "beta_raw_sesame.rds"))



# -----------------------------------------------------------------------------
# I.4 Per-sample QC and Supplemental Figure 2
# -----------------------------------------------------------------------------

# pOOBAH p-values from sesame are EPIC v2.0-aware

cat("Computing pOOBAH detection p-values...\n")
detP_list <- openSesame(idat_prefixes, prep = "Q",
                        func = pOOBAH, return.pval = TRUE,
                        platform = "EPICv2", manifest = epicv2_addr,
                        BPPARAM = SerialParam())

if (is.list(detP_list)) {
  common_p <- Reduce(intersect, lapply(detP_list, names))
  detP <- do.call(cbind, lapply(detP_list, function(v) v[common_p]))
} else {
  detP <- detP_list
}
colnames(detP) <- sample_sheet$Sample_Name

# Align rows: detP and beta_raw should have identical probe IDs
shared_probes <- intersect(rownames(beta_raw), rownames(detP))
beta_raw <- beta_raw[shared_probes, ]
detP     <- detP[shared_probes, ]

mean_detP <- colMeans(detP, na.rm = TRUE)
# (No standalone detection-p figure is written; panel A of Supplemental
#  Figure 2 carries this metric.)


# -----------------------------------------------------------------------------
# 4b. Extended per-sample QC metrics + supplementary figures
# -----------------------------------------------------------------------------
# Tier-1 metrics reviewers expect for EPIC arrays, computed on ALL samples
# (before any exclusion so the supplement documents who was dropped and why):
#   - bisulfite conversion efficiency (GCT score; ~1.0 complete, higher = worse)
#     with a minfi/wateRmelon::bscon fallback (% conversion from control probes)
#   - {M,U} signal intensity (proxy for DNA input / hybridization quality)
#   - genotype identity from rs (SNP) probes (sample-swap / replicate check)
# These are derived from the same QCDPB SigDFs (the identity check uses the rs
# probes carried in beta_raw). Mean detection p (section 4) remains the
# operative exclusion rule. We also retain an inferSex() prediction in the
# QC table: it is NOT used as a QC pass/fail metric (three subjects have no
# recorded sex, so concordance is uninformative) but as a way to impute the
# missing Sex covariate. Sex is not used in the adjusted CellDMC model in
# PART V; see the note there.
# NOTE: this reads the IDATs once more to obtain SigDF objects; for large
# cohorts you can refactor section 3 to keep the SigDF list and derive
# betas/detP from it in a single pass.

cat("\nComputing extended per-sample QC (GCT, intensity; inferSex for imputation)...\n")

# (i) SigDF objects (func = NULL returns the per-sample SigDF list)
sdf_list <- openSesame(idat_prefixes, prep = "QCDPB", func = NULL,
                       platform = "EPICv2", manifest = epicv2_addr,
                       BPPARAM = SerialParam())
if (!is.list(sdf_list)) sdf_list <- list(sdf_list)
names(sdf_list) <- sample_sheet$Sample_Name

# (ii) Predicted sex from betas (sesame uses its curated X/Y probe set).
#      Run on beta_raw (pre-filtering: sex-chromosome probes still present).
#      Retained for imputing the unknown-sex subjects, not for QC pass/fail.
pred_sex <- vapply(colnames(beta_raw), function(s) {
  b <- beta_raw[, s]; names(b) <- rownames(beta_raw)
  tryCatch(as.character(inferSex(b)), error = function(e) NA_character_)
}, character(1))

# (iii) GCT (bisulfite) and median {M,U} signal intensity per sample
# The tryCatch below used to swallow the GCT error silently, so a total
# failure of the metric was indistinguishable from every specimen passing.
# Collect the messages and report them.
gct_errs <- character(0)

qc_extra <- map_dfr(sample_sheet$Sample_Name, function(s) {
  sdf <- sdf_list[[s]]
  gct <- tryCatch(bisConversionControl(sdf),
                  error = function(e) {
                    gct_errs <<- c(gct_errs, conditionMessage(e))
                    NA_real_
                  })
  mu  <- tryCatch(signalMU(sdf),            error = function(e) NULL)
  tibble(
    Sample_Name    = s,
    gct_score      = gct,
    median_meth    = if (!is.null(mu)) median(mu$M, na.rm = TRUE) else NA_real_,
    median_unmeth  = if (!is.null(mu)) median(mu$U, na.rm = TRUE) else NA_real_
  )
})

if (length(gct_errs) > 0)
  cat("NOTE: sesame's GCT score (bisConversionControl) failed for",
      length(gct_errs), "of", nrow(sample_sheet), "specimens.\n",
      "      First error:", gct_errs[1], "\n",
      "      Falling back to wateRmelon::bscon below. If you would rather\n",
      "      report GCT, this is the error to resolve (it is usually a sesame\n",
      "      version that predates EPICv2 support).\n")

# (iv) Bisulfite conversion via minfi control probes (wateRmelon::bscon).
#       Reads the IDATs through minfi and summarizes the Illumina Bisulfite
#       Conversion I + II control probes into a per-sample conversion %.
#       This is the EPICv2 fallback for sesame's GCT score (which some sesame
#       versions cannot compute for EPICv2). Guarded: if minfi cannot read the
#       v2 IDATs or bscon does not support the array, bs_conversion stays NA.
bs_pct <- tryCatch({
  rg_bs <- minfi::read.metharray(idat_prefixes, force = TRUE)
  colnames(rg_bs) <- sample_sheet$Sample_Name
  setNames(as.numeric(wateRmelon::bscon(rg_bs)), sample_sheet$Sample_Name)
}, error = function(e) {
  cat("Bisulfite conversion (minfi/wateRmelon::bscon) unavailable:",
      conditionMessage(e), "\n")
  setNames(rep(NA_real_, nrow(sample_sheet)), sample_sheet$Sample_Name)
})


# (v) Genotype identity check from EPIC rs (SNP) probes
# -----------------------------------------------------------------------------
# EPIC v2.0 carries a panel of rs (SNP) probes whose beta reports genotype, not
# methylation: each sits near 0 (hom ref), 0.5 (het), or 1 (hom alt). Together
# they form a per-sample fingerprint, used two ways:
#   - swap / contamination: samples from DIFFERENT subjects that fingerprint-
#     match are a likely mislabel or cross-contamination.
#   - replicate confirmation: repeat / paired draws from the SAME subject should
#     match; a within-subject pair that does NOT match is flagged.
# Computed on beta_raw because the rs probes are still present here (section I.6
# removes them). The per-sample flag_identity is folded into qc_any_flag below;
# it is reported for investigation, not used as an automatic drop rule.

rs_idx   <- grep("^rs", rownames(beta_raw))
rs_betas <- beta_raw[rs_idx, , drop = FALSE]
rs_betas <- rs_betas[rowSums(is.na(rs_betas)) == 0, , drop = FALSE]  # complete probes only
cat("\nGenotype identity check: ", nrow(rs_betas), " complete rs probes across ",
    ncol(rs_betas), " samples\n", sep = "")

flag_identity <- setNames(rep(FALSE, ncol(beta_raw)), colnames(beta_raw))

if (nrow(rs_betas) >= RS_MIN_PROBES) {
  # Discretize betas to genotype calls {0,1,2}
  rs_geno <- matrix(1L, nrow(rs_betas), ncol(rs_betas), dimnames = dimnames(rs_betas))
  rs_geno[rs_betas < RS_HET_LOW]  <- 0L
  rs_geno[rs_betas > RS_HET_HIGH] <- 2L
  
  # Pairwise genotype concordance = fraction of matching calls
  samp <- colnames(rs_geno)
  concord <- matrix(NA_real_, length(samp), length(samp), dimnames = list(samp, samp))
  for (i in seq_along(samp))
    for (j in seq_along(samp))
      concord[i, j] <- mean(rs_geno[, i] == rs_geno[, j])
  
  subj_of <- setNames(sample_sheet$Subject_ID, sample_sheet$Sample_Name)
  
  # Classify every off-diagonal pair
  pair_tbl <- map_dfr(seq_len(length(samp) - 1L), function(i) {
    map_dfr((i + 1L):length(samp), function(j) {
      same <- unname(subj_of[samp[i]] == subj_of[samp[j]])
      cc   <- concord[i, j]
      verdict <- dplyr::case_when(
        !same & cc >= IDENT_MATCH_MIN ~ "UNEXPECTED_MATCH",     # possible swap/contamination
        same  & cc <  IDENT_MATCH_MIN ~ "UNEXPECTED_MISMATCH",  # possible mislabel
        TRUE                          ~ "ok")
      tibble(sample_i = samp[i], subject_i = unname(subj_of[samp[i]]),
             sample_j = samp[j], subject_j = unname(subj_of[samp[j]]),
             same_subject = same, concordance = cc, verdict = verdict)
    })
  })
  
  write_csv(pair_tbl, file.path(TBL_DIR, "genotype_identity_pairs.csv"))
  suspicious <- pair_tbl %>% filter(verdict != "ok")
  cat("  Suspicious identity pairs:", nrow(suspicious), "\n")
  if (nrow(suspicious) > 0) {
    print(suspicious %>% mutate(concordance = round(concordance, 3)))
    bad <- unique(c(suspicious$sample_i, suspicious$sample_j))
    flag_identity[bad] <- TRUE
  }
  
  # No heatmap is written. The concordance check is a safeguard against
  # sample swaps, not a reported result; every pair and its verdict is in
  # genotype_identity_pairs.csv, and any suspicious pair is printed above.
  cat("  Identity output: genotype_identity_pairs.csv\n")
} else {
  cat("  Skipped: fewer than", RS_MIN_PROBES,
      "complete rs probes available; cannot fingerprint reliably.\n")
}


# (vi) Assemble the per-sample QC table with flags
norm_sex <- function(x) {
  x <- toupper(trimws(as.character(x)))
  dplyr::case_when(x %in% c("F", "FEMALE") ~ "F",
                   x %in% c("M", "MALE")   ~ "M",
                   TRUE                    ~ NA_character_)
}

qc_df <- sample_sheet %>%
  dplyr::select(Sample_Name, Subject_ID, Sample_Group, Sex) %>%
  left_join(qc_extra, by = "Sample_Name") %>%
  mutate(
    mean_detP     = mean_detP[Sample_Name],
    predicted_sex_raw = pred_sex[Sample_Name],          # raw inferSex label
    predicted_sex = norm_sex(pred_sex[Sample_Name]),    # normalized to F/M
    recorded_sex  = norm_sex(Sex),                       # NA for unknown-sex subjects
    log2_meth     = log2(median_meth),
    log2_unmeth   = log2(median_unmeth),
    bs_conversion  = bs_pct[Sample_Name],
    flag_detp      = mean_detP > DETP_SAMPLE_MAX,
    flag_gct       = !is.na(gct_score) & gct_score > GCT_MAX,
    flag_intensity = ((log2_meth + log2_unmeth) / 2) < LOG2INT_MIN,
    flag_bscon     = !is.na(bs_conversion) & bs_conversion < BSCON_MIN,
    flag_identity  = flag_identity[Sample_Name],
    qc_any_flag    = flag_detp | flag_gct | flag_intensity | flag_bscon | flag_identity
  )

# Availability of each bisulfite metric. Defined here, before the flag summary
# is printed, because both flag_gct and flag_bscon are FALSE when their metric
# is NA: if NEITHER metric can be computed, bisulfite conversion is never
# assessed, yet qc_any_flag would still report every specimen as passing. That
# is a silent QC failure, so fail loudly instead.
has_gct   <- any(!is.na(qc_df$gct_score))
has_bscon <- any(!is.na(qc_df$bs_conversion))
if (!has_gct && !has_bscon)
  stop("Bisulfite conversion was not assessed: neither sesame's GCT score nor ",
       "wateRmelon::bscon could be computed. flag_gct and flag_bscon are ",
       "therefore FALSE for every specimen and qc_any_flag would pass the whole ",
       "cohort without testing conversion. Resolve one of the two before ",
       "reporting QC.")
cat("Bisulfite conversion assessed by:",
    if (has_gct) "sesame GCT score" else "wateRmelon::bscon (control probes)",
    "\n")

write_csv(qc_df, file.path(TBL_DIR, "sample_qc_metrics.csv"))
cat("Per-sample QC table written: sample_qc_metrics.csv\n")
cat("Samples flagged by >=1 QC metric:", sum(qc_df$qc_any_flag, na.rm = TRUE), "\n")
if (any(qc_df$qc_any_flag, na.rm = TRUE)) {
  print(qc_df %>% filter(qc_any_flag) %>%
          dplyr::select(Subject_ID, mean_detP, gct_score, bs_conversion,
                        log2_meth, log2_unmeth, flag_identity))
}


# -----------------------------------------------------------------------------
# 4c. Supplementary QC figure (publication-ready, 3-panel)
# -----------------------------------------------------------------------------
status_cols <- c(pass = "#56B4E9", flag = "#D55E00")

qc_plot_df <- qc_df %>%
  mutate(status = ifelse(qc_any_flag %in% TRUE, "flag", "pass"),
         label  = Subject_ID)

theme_qc <- theme_classic(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

# (A) Mean detection p per sample (the operative exclusion metric)
pDetP <- ggplot(qc_plot_df,
                aes(reorder(label, -mean_detP), mean_detP, fill = status)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = DETP_SAMPLE_MAX, linetype = "dashed", color = "grey30") +
  coord_flip() +
  scale_fill_manual(values = status_cols, name = NULL) +
  labs(title = "Mean Detection P (pOOBAH)",
       x = NULL, y = "Mean detection p") +
  theme_qc

# (B) Bisulfite conversion. Prefer sesame's GCT score; if unavailable for
#     EPICv2 (NA for all), use the minfi/wateRmelon conversion %; if both are
#     unavailable, emit an informative empty panel (detP is covered in panel A).
# has_gct / has_bscon are defined above, next to the flag logic they govern.
if (has_gct) {
  pBis <- ggplot(qc_plot_df,
                 aes(reorder(label, gct_score), gct_score, fill = status)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = GCT_MAX, linetype = "dashed", color = "grey30") +
    coord_flip() +
    scale_fill_manual(values = status_cols, name = NULL) +
    labs(title = "Bisulfite Conversion (GCT Score)",
         x = NULL, y = "GCT score (1.0 = complete conversion)") +
    theme_qc
} else if (has_bscon) {
  pBis <- ggplot(qc_plot_df,
                 aes(reorder(label, bs_conversion), bs_conversion, fill = status)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = BSCON_MIN, linetype = "dashed", color = "grey30") +
    coord_flip(ylim = c(min(70, min(qc_plot_df$bs_conversion, na.rm = TRUE)), 100)) +
    scale_fill_manual(values = status_cols, name = NULL) +
    labs(title = "Bisulfite Conversion (%, Control Probes)",
         x = NULL, y = "Conversion (%)") +
    theme_qc
} else {
  cat("NOTE: neither GCT nor wateRmelon bisulfite % was available;\n",
      "      panel B is left as a placeholder.\n")
  pBis <- ggplot() +
    annotate("text", x = 0, y = 0,
             label = "Bisulfite conversion\nmetric unavailable", size = 4) +
    theme_void() +
    labs(title = "Bisulfite Conversion (Unavailable)") +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
}

# (C) Signal intensity: log2 median methylated vs unmethylated
# The exclusion rule is mean(log2 M, log2 U) < LOG2INT_MIN, i.e. the line
# log2_meth = 2*LOG2INT_MIN - log2_unmeth. Every specimen in this dataset sits
# far above that line, so it fell outside the plotted range: panel C appeared
# to have no threshold while panels A and B both showed one. Draw the line only
# when it actually crosses the panel, and otherwise state the margin instead.
#
# The line spans y from 2L - max(x) to 2L - min(x) over the plotted x range, so
# it is on the panel only if that interval overlaps the plotted y range.
int_margin <- min((qc_plot_df$log2_meth + qc_plot_df$log2_unmeth) / 2,
                  na.rm = TRUE) - LOG2INT_MIN

# Compare against the PLOTTED range, not the data range: ggplot pads a
# continuous scale by 5% of the range at each end, so a line just outside the
# data can still cross the drawn panel.
.pad <- function(r) r + c(-1, 1) * 0.05 * diff(r)
.rx <- .pad(range(qc_plot_df$log2_unmeth, na.rm = TRUE))
.ry <- .pad(range(qc_plot_df$log2_meth,   na.rm = TRUE))
.int_line_lo <- 2 * LOG2INT_MIN - .rx[2]
.int_line_hi <- 2 * LOG2INT_MIN - .rx[1]
int_line_on_panel <- .int_line_lo <= .ry[2] && .int_line_hi >= .ry[1]

pInt <- ggplot(qc_plot_df, aes(log2_unmeth, log2_meth, color = status)) +
  geom_point(size = 2.4) +
  geom_text_repel(aes(label = label), size = 2.6, show.legend = FALSE,
                  min.segment.length = 0,        # draw a stem for EVERY label
                  segment.size = 0.3, segment.color = "grey55",
                  box.padding = 0.5, point.padding = 0.25,
                  max.overlaps = Inf, seed = 1) +
  scale_color_manual(values = status_cols, name = NULL) +
  labs(title = "Signal Intensity",
       x = expression(log[2]~median~unmethylated),
       y = expression(log[2]~median~methylated),
       # Only claim a clean pass when nothing is actually flagged: the line is
       # also off-panel when every specimen falls BELOW it, and the margin
       # would then print negative.
       subtitle = if (int_line_on_panel ||
                      any(qc_plot_df$flag_intensity, na.rm = TRUE)) NULL else
         sprintf("All specimens exceed the mean log2 intensity floor of %.1f by at least %.2f",
                 LOG2INT_MIN, int_margin)) +
  theme_qc +
  theme(plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5))

if (int_line_on_panel)
  pInt <- pInt + geom_abline(slope = -1, intercept = 2 * LOG2INT_MIN,
                             linetype = "dashed", color = "grey30")

# Sex concordance is intentionally not shown: three subjects have no recorded
# sex, so a predicted-vs-recorded concordance check is uninformative. The
# inferSex() prediction is retained in the QC table for covariate imputation.

# The subtitle must name the metric panel B actually shows. GCT is NA for all
# specimens on EPICv2 in this dataset, so panel B falls through to the
# wateRmelon conversion percentage; the previous subtitle advertised a GCT
# threshold that was never applied and that no panel displayed.
# Both rules are in force whenever their metric is non-NA (qc_any_flag ORs
# flag_gct and flag_bscon); has_gct/has_bscon only decide which one panel B
# PLOTS. List every threshold that can actually flag a specimen.
bis_thresh_txt <- paste(c(
  if (has_gct)   paste0("GCT \u2264 ", GCT_MAX),
  if (has_bscon) paste0("bisulfite conversion \u2265 ", BSCON_MIN, "%"),
  if (!has_gct && !has_bscon) "bisulfite conversion not assessed"),
  collapse = ", ")
cat("Supplemental Figure 2 panel B metric:",
    if (has_gct) "GCT score" else if (has_bscon) "conversion %" else "none",
    "\n")

supp_qc <- (pDetP | pBis) / pInt +
  plot_annotation(
    tag_levels = "A",
    title = "Supplemental Figure 2. Per-Sample Array Quality Control (EPIC v2.0)",
    subtitle = paste0("Flag thresholds: mean detP \u2264 ", DETP_SAMPLE_MAX,
                      ", ", bis_thresh_txt,
                      ", mean log2 intensity \u2265 ", LOG2INT_MIN),
    theme = theme(plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
                  plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0.5))) &
  theme(plot.tag = element_text(face = "bold", size = 13))

ggsave(file.path(FIG_DIR, "SuppFig2_QC_panel.pdf"), supp_qc,
       width = 11, height = 8, device = cairo_pdf, dpi = 600)
ggsave(file.path(FIG_DIR, "SuppFig2_QC_panel.png"), supp_qc,
       width = 11, height = 8, dpi = 600)
cat("Supplemental Figure 2 written: SuppFig2_QC_panel.pdf / .png\n")


bad_samples <- names(mean_detP)[mean_detP > 0.05]
if (length(bad_samples) > 0) {
  cat("Dropping", length(bad_samples), "samples with mean detP > 0.05:",
      paste(bad_samples, collapse = ", "), "\n")
  keep <- !(colnames(beta_raw) %in% bad_samples)
  beta_raw <- beta_raw[, keep]
  detP     <- detP[, keep]
  sample_sheet <- sample_sheet %>% filter(!(Sample_Name %in% bad_samples))
}

cat("Samples retained after QC:", ncol(beta_raw), "\n")


# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# I.5 M-values
# -----------------------------------------------------------------------------

# Use small offset so beta values at the boundary don't blow up.

eps <- 1e-6
beta_clipped <- pmin(pmax(beta_raw, eps), 1 - eps)
m_raw <- log2(beta_clipped / (1 - beta_clipped))


# -----------------------------------------------------------------------------
# I.6 Probe filtering
# -----------------------------------------------------------------------------
# (a) Detection p > 0.01 in any sample
# (b) Cross-reactive / design-flagged probes via the Zhou KYCG EPICv2 mask
# (c) Sex chromosomes
# (d) Invariant probes: beta > 0.8 in all samples or beta < 0.2 in all samples
#
# Every count is captured in `probe_counts` as it is computed, and the
# accounting table written in section I.9 is built from that list rather than
# transcribed from a log file. The detection-p and variance filters are
# independent pass counts against the full array, not a sequential cascade;
# the retained set is their intersection.

probe_counts <- list(array = nrow(beta_raw))

# (a)
keep_detP <- rowSums(detP > 0.01, na.rm = TRUE) == 0
probe_counts$pass_detp <- sum(keep_detP)
cat("Passing detP filter:", sum(keep_detP), "/", length(keep_detP), "\n")

# (b) Cross-reactive / design-flagged probes from the KYCG EPICv2 mask.
#     This sesameDataGet resolves via the alternate host (SESAMEDATA_USE_ALT,
#     set in section I.1), returning the curated mask rather than the cold
#     ExperimentHub copy. The object is a named list of mask sub-DBs (SNP,
#     mapping, non-unique, and so on); unlist + unique collapses them to the
#     full set of design-flagged probe IDs.
epicv2_mask <- sesameDataGet("KYCG.EPICv2.Mask.20230314")
cross_reactive <- unique(unlist(epicv2_mask, use.names = FALSE))
probe_counts$cross_reactive <- sum(rownames(beta_raw) %in% cross_reactive)
cat("Cross-reactive probes (KYCG EPICv2 mask):", length(cross_reactive),
    "of which", probe_counts$cross_reactive, "are on this array\n")

# (c) Sex chromosomes
annot_epic <- as.data.frame(
  getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
) %>%
  rownames_to_column("probeID")
sex_chr_probes <- annot_epic$probeID[annot_epic$chr %in% c("chrX", "chrY")]
probe_counts$sex_chr <- sum(rownames(beta_raw) %in% sex_chr_probes)
cat("Sex-chromosome probes on this array:", probe_counts$sex_chr, "\n")

# (d) Invariant
all_high <- rowSums(beta_raw > 0.8, na.rm = TRUE) == ncol(beta_raw)
all_low  <- rowSums(beta_raw < 0.2, na.rm = TRUE) == ncol(beta_raw)
keep_var <- !(all_high | all_low)
probe_counts$pass_varsat <- sum(keep_var)
cat("Passing variance/saturation filter:", sum(keep_var), "/", length(keep_var), "\n")

keep_probes <- keep_detP & keep_var &
  !(rownames(beta_raw) %in% cross_reactive) &
  !(rownames(beta_raw) %in% sex_chr_probes)
probe_counts$after_all <- sum(keep_probes)
cat("Probes retained after all filters:", sum(keep_probes), "\n")

beta_filt <- beta_raw[keep_probes, ]
m_filt    <- m_raw[keep_probes, ]


# -----------------------------------------------------------------------------
# I.7 Collapse EPIC v2.0 replicate probes to base IDs
# -----------------------------------------------------------------------------

# collapse replicates by mean before deconvolution. avereps() is the limma
# idiom for this.

cat("Collapsing EPIC v2 replicate probes to base IDs...\n")
base_ids <- sub("_[A-Z0-9]+$", "", rownames(beta_filt))
beta_collapsed <- limma::avereps(beta_filt, ID = base_ids)
m_collapsed    <- limma::avereps(m_filt,    ID = base_ids)

cat("After collapsing replicates:\n")
cat("  Unique base probes:", nrow(beta_collapsed), "\n")
cat("  Reduction from:    ", nrow(beta_filt), "\n")

# Drop any probes with NA or non-finite values in any sample. prcomp() and
# limma cannot tolerate these; a small number of probes (~0.3%) typically
# pick up NAs from the per-sample detection-P filter without being entirely
# removed.
finite_probes <- rowSums(is.na(m_collapsed) | !is.finite(m_collapsed)) == 0
n_dropped <- nrow(m_collapsed) - sum(finite_probes)
if (n_dropped > 0) {
  cat("Dropping", n_dropped, "probes with NA / non-finite values from",
      "m_collapsed and beta_collapsed\n")
  m_collapsed    <- m_collapsed[finite_probes, ]
  beta_collapsed <- beta_collapsed[finite_probes, ]
}
cat("Final analytic matrix:", nrow(m_collapsed), "probes x",
    ncol(m_collapsed), "samples\n")

saveRDS(beta_collapsed, file.path(RDS_DIR, "beta_filtered_collapsed.rds"))
saveRDS(m_collapsed,    file.path(RDS_DIR, "m_filtered_collapsed.rds"))



probe_counts$after_collapse <- nrow(beta_collapsed) + n_dropped
probe_counts$dropped_nonfinite_collapse <- n_dropped


# -----------------------------------------------------------------------------
# I.8 Cell-type deconvolution, full cohort (estimation)
# -----------------------------------------------------------------------------

# ESTIMATION ONLY. Every result derived from these estimates is reported in
# PART III, in manuscript order. The step runs here because the final probe set
# depends on the composition adjustment below, and the adjustment depends on
# these fractions. See the ordering note in the file header.

#   ref1 = centEpiFibIC (Epithelial / Fibroblast / Immune Cell)
#   ref2 = centDHSbloodDMC (B, NK, CD4T, CD8T, Mono, Neutro, Eosino)
#   h.CT.idx = 3 = the IC column of ref1 is subdivided by ref2

data(centEpiFibIC.m,    package = "EpiDISH", envir = environment())
data(centDHSbloodDMC.m, package = "EpiDISH", envir = environment())

cat("Overlap with EpiDISH references:\n")
cat("  Ref1 (Epi/Fib/IC):  ",
    sum(rownames(beta_collapsed) %in% rownames(centEpiFibIC.m)),
    "of", nrow(centEpiFibIC.m), "\n")
cat("  Ref2 (immune sub): ",
    sum(rownames(beta_collapsed) %in% rownames(centDHSbloodDMC.m)),
    "of", nrow(centDHSbloodDMC.m), "\n")

cat("\nRunning EpiDISH hepidish() with RPC...\n")
hepidish_res <- hepidish(
  beta.m   = beta_collapsed,
  ref1.m   = centEpiFibIC.m,
  ref2.m   = centDHSbloodDMC.m,
  h.CT.idx = 3,
  method   = "RPC"
)

cell_props <- as.data.frame(hepidish_res) %>%
  rownames_to_column("Sample_Name")

# Relabel to Subject_ID for the human-readable output table only. Keep the
# working `cell_props` keyed by Sample_Name, since the adjustment below and
# every cohort analysis match on colnames(m_collapsed), which are Sample_Names.
cell_props_out <- cell_props %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID),
            by = "Sample_Name") %>%
  dplyr::select(Subject_ID, everything(), -Sample_Name)

cat("\nDeconvolution results (rounded to 3dp):\n")
print(cell_props_out %>% mutate(across(where(is.numeric), ~round(.x, 3))))
write_csv(cell_props_out, file.path(TBL_DIR, "cell_composition_estimates.csv"))

# Cell-type column names. Every downstream section indexes the fraction
# matrices with this vector, so it is defined once, here, from the hepidish
# output rather than being re-derived per section.
cell_cols <- setdiff(names(cell_props), "Sample_Name")
cat("Cell types estimated:", paste(cell_cols, collapse = ", "), "\n")


# -----------------------------------------------------------------------------
# I.9 Composition adjustment, full cohort, and the final analytic matrix
# -----------------------------------------------------------------------------

# Use limma::removeBatchEffect() with cell proportions as continuous
# covariates. Three preparation steps avoid the NA / rank-deficiency problems
# that bite removeBatchEffect when there are few samples and collinear
# covariates:
#   (a) Drop the largest-mean cell-type column to avoid the sum-to-1
#       collinearity trap.
#   (b) Drop near-zero-variance cell-type columns (e.g. Eosino and CD8T which
#       are zero in most nasal samples) since they contribute no useful
#       adjustment and just inflate the design rank deficiency.
#   (c) After removeBatchEffect, strip any probes that picked up NAs or
#       non-finite values from m_collapsed, beta_collapsed, AND m_adjusted in
#       lockstep so every downstream matrix stays aligned and clean.

prop_mat <- as.matrix(
  cell_props[match(colnames(m_collapsed), cell_props$Sample_Name),
             cell_cols, drop = FALSE]
)

# (a) reference column to drop (largest mean)
drop_col <- cell_cols[which.max(colMeans(prop_mat))]
cat("Dropping reference cell type (largest mean):", drop_col,
    "(", round(colMeans(prop_mat)[drop_col], 3), ")\n")

# (b) drop near-zero-variance columns. VAR_THRESHOLD (1e-4, stated in the
#     Methods) catches CD8T and Eosino when they are zero or near-zero in all
#     samples. The same constant is used in every cohort and every strategy.
candidate_cols <- setdiff(cell_cols, drop_col)
col_vars <- apply(prop_mat[, candidate_cols, drop = FALSE], 2, var)
near_zero_var_cols <- names(col_vars)[col_vars < VAR_THRESHOLD]
if (length(near_zero_var_cols) > 0) {
  cat("Dropping near-zero-variance composition columns:",
      paste(near_zero_var_cols, collapse = ", "), "\n")
}
keep_cols <- setdiff(candidate_cols, near_zero_var_cols)
prop_for_adj <- prop_mat[, keep_cols, drop = FALSE]
cat("Composition covariates retained for adjustment:",
    paste(keep_cols, collapse = ", "),
    "(", ncol(prop_for_adj), "columns)\n")

# Run the adjustment
m_adjusted <- removeBatchEffect(m_collapsed, covariates = prop_for_adj)

# (c) Post-adjustment NA cleanup. removeBatchEffect can introduce NAs at
#     probes whose residuals cannot be uniquely determined given the
#     covariate structure. We strip those probes from all three matrices so
#     downstream PCA, clustering, and limma never see them.
ok_after_adj <- rowSums(is.na(m_adjusted) | !is.finite(m_adjusted)) == 0
n_dropped_adj <- nrow(m_adjusted) - sum(ok_after_adj)
if (n_dropped_adj > 0) {
  cat("Dropping", n_dropped_adj, "probes with NAs introduced by",
      "removeBatchEffect from m_collapsed, beta_collapsed, and m_adjusted\n")
  m_adjusted     <- m_adjusted[ok_after_adj, ]
  m_collapsed    <- m_collapsed[ok_after_adj, ]
  beta_collapsed <- beta_collapsed[ok_after_adj, ]
}
cat("Final analytic matrices (post-adjustment cleanup):", nrow(m_adjusted),
    "probes x", ncol(m_adjusted), "samples\n")


saveRDS(m_adjusted, file.path(RDS_DIR, "m_adjusted_composition.rds"))

# -----------------------------------------------------------------------------
# I.10 Probe filtering accounting (Supplemental Table 2)
# -----------------------------------------------------------------------------
# Reported in the first paragraph of the Results. Built from the counts
# recorded during filtering rather than transcribed, so the table cannot
# disagree with the matrix actually analyzed.

probe_counts$dropped_nonfinite_adj <- n_dropped_adj
probe_counts$final <- nrow(m_collapsed)

stopifnot(probe_counts$after_collapse - probe_counts$dropped_nonfinite_collapse -
            probe_counts$dropped_nonfinite_adj == probe_counts$final)

pct_of_array <- function(n) sprintf("%.1f", 100 * n / probe_counts$array)

probe_filtering <- tibble::tribble(
  ~step, ~basis, ~probes,
  "Probes on the EPIC v2 array",
  "Starting set",                        probe_counts$array,
  "Passing detection p-value filter",
  "Independent filter, of starting set", probe_counts$pass_detp,
  "Masked as cross-reactive (KYCG EPIC v2)",
  "Probes removed",                      probe_counts$cross_reactive,
  "On the X or Y chromosome",
  "Probes removed",                      probe_counts$sex_chr,
  "Passing variance and saturation filter",
  "Independent filter, of starting set", probe_counts$pass_varsat,
  "Retained after all filters",
  "Intersection of the filters above",   probe_counts$after_all,
  "After collapsing replicate probes to base identifiers",
  "Cumulative",                          probe_counts$after_collapse,
  "Removed for missing or non-finite values after collapsing",
  "Probes removed",                      probe_counts$dropped_nonfinite_collapse,
  "Removed for non-finite residuals after composition adjustment",
  "Probes removed",                      probe_counts$dropped_nonfinite_adj,
  "Final analytic matrix",
  "Cumulative",                          probe_counts$final
) %>%
  mutate(pct_of_array = pct_of_array(probes))

cat("\nProbe filtering accounting:\n")
print(as.data.frame(probe_filtering))

save_pub_table(
  probe_filtering %>%
    transmute(`Filtering step` = step,
              `Basis`          = basis,
              `Probes`         = format(probes, big.mark = ","),
              `% of array`     = pct_of_array),
  file.path(TBL_DIR, "SuppTable2_probe_filtering"),
  caption = paste0("Supplemental Table 2. Probe filtering accounting for the ",
                   "Infinium MethylationEPIC v2.0 array."),
  note = paste0("The detection p-value and variance/saturation filters are ",
                "reported as independent pass counts against the full array ",
                "rather than as a sequential cascade; the retained count of ",
                format(probe_counts$after_all, big.mark = ","), " probes is ",
                "the intersection of these filters together with removal of ",
                "cross-reactive and sex-chromosome probes. EPIC v2.0 includes ",
                "replicate probes targeting the same CpG, which were ",
                "collapsed to unique base identifiers. ",
                format(probe_counts$dropped_nonfinite_collapse, big.mark = ","),
                " probes carried a missing or non-finite value after ",
                "collapsing and were removed from all matrices so that ",
                "downstream analyses remained aligned. The same rule was ",
                "applied after composition adjustment and removed ",
                ifelse(probe_counts$dropped_nonfinite_adj == 0, "no further probes",
                       paste(format(probe_counts$dropped_nonfinite_adj,
                                    big.mark = ","), "further probes")),
                ", leaving a final analytic matrix of ",
                format(probe_counts$final, big.mark = ","),
                " probes across ", ncol(m_collapsed), " specimens."))


# -----------------------------------------------------------------------------
# I.11 The three analytic cohorts
# -----------------------------------------------------------------------------
# Full cohort      all 24 specimens, PARDS and control (primary analysis)
# PARDS only       the 21 PARDS specimens, controls removed (Reviewer 1
#                  comment 1; Reviewer 3 comment 3)
# New PARDS only   the 13 PARDS specimens not analyzed in Williams et al.,
#                  Respir Res 2022;23:181 (Reviewer 3 comment 4)
#
# Every downstream analysis is re-derived within each cohort: deconvolution,
# the composition adjustment, the top 5% variance filter, PCA, and clustering.
# Only the probe set is shared, because probe filtering was performed at the
# cohort level.

stopifnot(all(c(CASE_LABEL, CONTROL_LABEL) %in% sample_sheet$Sample_Group))
stopifnot(all(DUPLICATE_SUBJECTS %in% sample_sheet$Subject_ID))

pards_samples <- sample_sheet$Sample_Name[sample_sheet$Sample_Group == CASE_LABEL]
pards_samples <- intersect(colnames(m_collapsed), pards_samples)
n_pards       <- length(pards_samples)

new_ss <- sample_sheet %>%
  filter(Sample_Name %in% colnames(beta_collapsed),
         !Subject_ID %in% DUPLICATE_SUBJECTS)
new_pards <- new_ss$Sample_Name[new_ss$Sample_Group == CASE_LABEL]
n_new     <- length(new_pards)

cat("\nAnalytic cohorts\n")
cat("  Full cohort:    ", ncol(m_collapsed), "specimens (",
    sum(sample_sheet$Sample_Group == CASE_LABEL), "PARDS,",
    sum(sample_sheet$Sample_Group == CONTROL_LABEL), "control )\n")
cat("  PARDS only:     ", n_pards, "specimens\n")
cat("  New PARDS only: ", n_new, "specimens (",
    length(DUPLICATE_SUBJECTS), "reused specimens excluded, of which",
    sum(grepl("^Control", DUPLICATE_SUBJECTS)), "are controls )\n")
cat("  Controls in the new-specimen subset:",
    sum(new_ss$Sample_Group == CONTROL_LABEL),
    "(a single control cannot inform a clustering, so the new-specimen\n",
    "  analyses are restricted to PARDS)\n")
stopifnot(n_pards >= 4, n_new >= 6)

beta_pards        <- beta_collapsed[, pards_samples, drop = FALSE]
m_collapsed_pards <- m_collapsed[,    pards_samples, drop = FALSE]
beta_new          <- beta_collapsed[, new_pards, drop = FALSE]
m_new             <- m_collapsed[,    new_pards, drop = FALSE]

COHORT_LABELS <- c(paste0("Full cohort (n = ", ncol(m_collapsed), ")"),
                   paste0("PARDS only (n = ", n_pards, ")"),
                   paste0("New PARDS only (n = ", n_new, ")"))

COHORTS <- list(
  list(label = COHORT_LABELS[1], beta = beta_collapsed,    m = m_collapsed),
  list(label = COHORT_LABELS[2], beta = beta_pards,        m = m_collapsed_pards),
  list(label = COHORT_LABELS[3], beta = beta_new,          m = m_new))



# -----------------------------------------------------------------------------
# I.12 Table 1: cohort characteristics
# -----------------------------------------------------------------------------

# Reviewer 2 objects to Mann-Whitney p-values when one group has n = 3, and
# Reviewer 1 asks that all group comparisons be removed. This table is
# therefore purely descriptive: no hypothesis tests, no p-values, no
# significance language. Continuous variables are median [IQR]; categorical
# variables are n (%) within the column.
#
# Each entry in VAR_MAP lists candidate column names; the first one present in
# sample_sheet is used, so the block tolerates naming differences without
# editing. Anything not found is skipped with a warning.

cat("\n=== PART I.12. Descriptive cohort table (Table 1) ===\n")
cat("Columns available in sample_sheet:\n")
print(names(sample_sheet))

# Race is collapsed to White vs Non-White. Unknown is retained as its own row
# rather than folded into either category or dropped silently.
recode_race <- function(x) {
  x <- trimws(as.character(x))
  out <- dplyr::case_when(
    is.na(x)                              ~ NA_character_,
    grepl("^unk|^unknown|^not report", x, ignore.case = TRUE) ~ "Unknown",
    grepl("^white|^caucasian", x, ignore.case = TRUE)         ~ "White",
    TRUE                                                       ~ "Non-White")
  factor(out, levels = c("White", "Non-White", "Unknown"))
}

# PARDS severity is ordered so the table reads mild to severe rather than
# alphabetically. Controls have no PARDS category by definition.
recode_severity <- function(x) {
  x <- trimws(as.character(x))
  out <- dplyr::case_when(
    is.na(x)                                  ~ NA_character_,
    grepl("^none|^no ", x, ignore.case = TRUE) ~ "None",
    grepl("^mild", x, ignore.case = TRUE)      ~ "Mild",
    grepl("^mod",  x, ignore.case = TRUE)      ~ "Moderate",
    grepl("^sev",  x, ignore.case = TRUE)      ~ "Severe",
    TRUE                                       ~ x)
  factor(out, levels = c("None", "Mild", "Moderate", "Severe"))
}

VAR_MAP <- list(
  list(label = "Age, years", type = "cont",
       cols = c("Age", "Age_years", "age")),
  list(label = "Sex", type = "cat",
       cols = c("Sex", "sex", "Gender")),
  list(label = "Race", type = "cat",
       cols = c("Race", "race"), recode = recode_race),
  list(label = "Principal comorbidity", type = "cat",
       cols = c("Principal_Comorbidity", "Principal_Comorbidities", "Comorbidity")),
  list(label = "Acute condition", type = "cat",
       cols = c("Acute_Condition", "Acute_Dx", "Acute_Diagnosis", "Admission_Diagnosis")),
  list(label = "PARDS severity", type = "cat",
       cols = c("PARDS_Severity", "PARDS_Category", "Highest_PARDS_Category",
                "PARDS_severity", "Severity", "OI_Category"),
       recode = recode_severity),
  list(label = "PELOD-2 score", type = "cont",
       cols = c("PELOD2", "PELOD_2", "PELOD", "Highest_PELOD2", "PELOD2_max")),
  list(label = "Ventilator-free days", type = "cont",
       cols = c("VFD", "Ventilator_Free_Days", "VFD_28", "Vent_Free_Days")),
  list(label = "Outcome", type = "cat",
       cols = c("Outcome", "Mortality", "Survival", "Vital_Status"))
)

pick_col <- function(cands) {
  hit <- cands[cands %in% names(sample_sheet)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

resolved <- vapply(VAR_MAP, function(v) pick_col(v$cols), character(1))
names(resolved) <- vapply(VAR_MAP, `[[`, character(1), "label")
cat("\nResolved column names:\n"); print(resolved)
if (any(is.na(resolved)))
  cat("\nNOT FOUND (rows skipped):",
      paste(names(resolved)[is.na(resolved)], collapse = ", "),
      "\nAdd the correct name to the cols vector in VAR_MAP.\n")

tbl_dat <- sample_sheet %>%
  filter(Sample_Name %in% colnames(beta_collapsed)) %>%
  mutate(Group = factor(if_else(Sample_Group == CASE_LABEL, "PARDS", "Control"),
                        levels = c("PARDS", "Control")))

n_by_group <- table(tbl_dat$Group)
cat("\nSpecimens in table: PARDS =", n_by_group[["PARDS"]],
    ", Control =", n_by_group[["Control"]], "\n")

fmt_cont <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("\u2014")
  sprintf("%.1f [%.1f, %.1f]", median(x), quantile(x, 0.25), quantile(x, 0.75))
}
fmt_cat <- function(x, level) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("\u2014")
  n <- sum(as.character(x) == level)
  sprintf("%d (%.0f)", n, 100 * n / length(x))
}

groups <- levels(tbl_dat$Group)

build_rows <- function(spec) {
  col <- pick_col(spec$cols)
  if (is.na(col)) return(NULL)
  v <- tbl_dat[[col]]
  if (!is.null(spec$recode)) v <- spec$recode(v)

  if (spec$type == "cont") {
    v <- suppressWarnings(as.numeric(v))
    row <- tibble(Characteristic = paste0(spec$label, ", median [IQR]"),
                  .row_type = "var")
    for (g in groups) row[[g]] <- fmt_cont(v[tbl_dat$Group == g])
    n_miss <- sum(is.na(v))
    if (n_miss > 0)
      row$Characteristic <- paste0(row$Characteristic, " (missing n = ", n_miss, ")")
    return(dplyr::select(row, Characteristic, all_of(groups), .row_type))
  }

  lv <- if (is.factor(v)) levels(v)[levels(v) %in% as.character(v[!is.na(v)])]
        else sort(unique(as.character(v[!is.na(v)])))

  header <- tibble(Characteristic = paste0(spec$label, ", n (%)"), .row_type = "var")
  for (g in groups) header[[g]] <- ""
  header <- dplyr::select(header, Characteristic, all_of(groups), .row_type)

  body <- map_dfr(lv, function(l) {
    r <- tibble(Characteristic = l, .row_type = "level")
    for (g in groups) r[[g]] <- fmt_cat(v[tbl_dat$Group == g], l)
    dplyr::select(r, Characteristic, all_of(groups), .row_type)
  })

  if (sum(is.na(v)) > 0) {
    r <- tibble(Characteristic = "Missing", .row_type = "level")
    for (g in groups) r[[g]] <- as.character(sum(is.na(v[tbl_dat$Group == g])))
    body <- bind_rows(body, dplyr::select(r, Characteristic, all_of(groups), .row_type))
  }
  bind_rows(header, body)
}

table1_long <- map_dfr(VAR_MAP, build_rows)

row_type <- table1_long$.row_type
demographic_table <- table1_long %>% dplyr::select(-.row_type)

names(demographic_table)[names(demographic_table) == "PARDS"] <-
  sprintf("PARDS (n = %d)", n_by_group[["PARDS"]])
names(demographic_table)[names(demographic_table) == "Control"] <-
  sprintf("Control (n = %d)", n_by_group[["Control"]])

cat("\nDescriptive cohort table:\n")
print(as.data.frame(demographic_table), row.names = FALSE)

write_csv(demographic_table, file.path(TBL_DIR, "Table1_cohort_descriptive.csv"))

# -----------------------------------------------------------------------------
# Rendering: variable names bold and flush left, levels indented, and a rule
# above each variable block so categories do not run together.
# -----------------------------------------------------------------------------
tbl_note <- paste0(
  "Values are median [interquartile range] for continuous variables and ",
  "n (%) for categorical variables, with percentages calculated within each ",
  "column among non-missing observations. Race was collapsed to White versus ",
  "Non-White. PARDS severity applies only to PARDS subjects. Ventilator-free ",
  "days (VFD) was calculated as VFD = 28 minus the number of days intubated. ",
  "No statistical comparisons were performed between groups. Given that the ",
  "control group comprises ", n_by_group[["Control"]],
  " specimens, between-group differences should be regarded as descriptive ",
  "only and no inference is intended.")

if (requireNamespace("flextable", quietly = TRUE) &&
    requireNamespace("officer", quietly = TRUE)) {

  hdr_i <- which(row_type == "var")
  lvl_i <- which(row_type == "level")
  rule  <- officer::fp_border(color = "black", width = 0.75)
  hair  <- officer::fp_border(color = "grey60", width = 0.5)

  ft <- flextable::flextable(demographic_table)
  ft <- flextable::theme_booktabs(ft)
  ft <- flextable::font(ft, fontname = "Arial", part = "all")
  ft <- flextable::fontsize(ft, size = 9, part = "all")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::align(ft, j = 1, align = "left", part = "all")
  ft <- flextable::align(ft, j = 2:ncol(demographic_table), align = "center", part = "all")
  # variable rows: bold, flush left, rule above each block after the first
  ft <- flextable::bold(ft, i = hdr_i, j = 1, part = "body")
  if (length(hdr_i) > 1)
    ft <- flextable::hline(ft, i = hdr_i[-1] - 1, border = hair, part = "body")
  # level rows: indented so they read as belonging to the block above
  ft <- flextable::padding(ft, i = lvl_i, j = 1, padding.left = 18, part = "body")
  ft <- flextable::padding(ft, padding.top = 2, padding.bottom = 2, part = "body")
  ft <- flextable::hline_top(ft, border = rule, part = "header")
  ft <- flextable::hline_bottom(ft, border = rule, part = "body")
  ft <- flextable::autofit(ft)

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Table 1. Characteristics of the analyzed cohort.",
                               style = "heading 2")
  doc <- flextable::body_add_flextable(doc, ft)
  doc <- officer::body_add_par(doc, tbl_note, style = "Normal")
  print(doc, target = file.path(TBL_DIR, "Table1_cohort_descriptive.docx"))
  cat("  Wrote publication table: Table1_cohort_descriptive.docx\n")
} else {
  cat("  (flextable/officer not installed -> CSV only)\n")
}

saveRDS(list(table = demographic_table, row_type = row_type, note = tbl_note),
        file.path(RDS_DIR, "checkpoint_table1_descriptive.rds"))
cat("\nPART I complete:", format(Sys.time()), "\n")


# =============================================================================
# PART II. UNADJUSTED PCA IN ALL THREE COHORTS  (Results: Dimensional Reduction)
# =============================================================================


# The top 5% most variable CpGs are reselected within each cohort. The variance
# explained is reported here and again, beside the post-adjustment values, in
# Supplemental Table 3 (PART IV). PC1 is interrogated against cell composition
# in PART III.

unadj_pca <- lapply(COHORTS, function(s) pca_top5(s$m))
names(unadj_pca) <- COHORT_LABELS

unadj_variance <- map_dfr(COHORT_LABELS, function(lab) {
  ve <- unadj_pca[[lab]]$var_exp
  tibble(cohort      = lab,
         n           = ncol(unadj_pca[[lab]]$m_top),
         probes_used = nrow(unadj_pca[[lab]]$m_top),
         pc1_pct     = ve[1] * 100,
         pc2_pct     = ve[2] * 100,
         pc3_pct     = if (length(ve) >= 3) ve[3] * 100 else NA_real_)
})

cat("\nUnadjusted PCA, variance explained by the leading components:\n")
print(as.data.frame(unadj_variance %>% mutate(across(where(is.numeric), ~round(.x, 1)))))
write_csv(unadj_variance, file.path(TBL_DIR, "unadjusted_pca_variance.csv"))

cat("\n--- For the Results text ---\n")
cat(sprintf(paste0("PC1 explained %.1f%%, %.1f%%, and %.1f%% of the variance ",
                   "in the full, PARDS only, and new PARDS only cohorts; PC2 ",
                   "explained %.1f%%, %.1f%%, and %.1f%%.\n"),
            unadj_variance$pc1_pct[1], unadj_variance$pc1_pct[2],
            unadj_variance$pc1_pct[3], unadj_variance$pc2_pct[1],
            unadj_variance$pc2_pct[2], unadj_variance$pc2_pct[3]))

# ---- Figure 1: unadjusted PCA, full cohort -----------------------------------
var_exp_raw <- unadj_pca[[1]]$var_exp
pca_raw     <- unadj_pca[[1]]$pca

fig1_df <- as.data.frame(pca_raw$x[, 1:2]) %>%
  rownames_to_column("Sample_Name") %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID, Sample_Group),
            by = "Sample_Name") %>%
  mutate(label = fig_label(Subject_ID))

fig1_pca <- ggplot(fig1_df, aes(PC1, PC2, color = Sample_Group, label = label)) +
  geom_point(size = 3.4) +
  geom_text_repel(size = 2.7, show.legend = FALSE,
                  min.segment.length = 0,
                  segment.size = 0.3, segment.color = "grey55",
                  box.padding = 0.5, point.padding = 0.25,
                  max.overlaps = Inf, seed = 1) +
  scale_color_manual(values = c(ARDS = FIG1_ORANGE, Control = FIG1_BLUE),
                     name = "Group") +
  labs(x = paste0("PC1 (", round(var_exp_raw[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_exp_raw[2] * 100, 1), "%)"),
       title = paste0("Figure 1. Unadjusted Methylome PCA, Full Cohort (n = ",
                      nrow(fig1_df), ")")) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(FIG_DIR, "Figure1_PCA_unadjusted.pdf"), fig1_pca,
       width = 7, height = 5, device = cairo_pdf, dpi = 600)
ggsave(file.path(FIG_DIR, "Figure1_PCA_unadjusted.png"), fig1_pca,
       width = 7, height = 5, dpi = 600)
cat("Figure 1 written: Figure1_PCA_unadjusted.pdf / .png\n")



# =============================================================================
# PART III. CELL-TYPE DECONVOLUTION IN ALL THREE COHORTS
# =============================================================================


# -----------------------------------------------------------------------------
# III.1 Deconvolution within the PARDS-only cohort
# -----------------------------------------------------------------------------

# Same two-stage reference and method as master section 8:
#   ref1 = centEpiFibIC (Epithelial / Fibroblast / Immune Cell)
#   ref2 = centDHSbloodDMC (B, NK, CD4T, CD8T, Mono, Neutro, Eosino)
#   h.CT.idx = 3 = the IC column of ref1 is subdivided by ref2

data(centEpiFibIC.m,    package = "EpiDISH", envir = environment())
data(centDHSbloodDMC.m, package = "EpiDISH", envir = environment())

n_ref1_overlap <- sum(rownames(beta_pards) %in% rownames(centEpiFibIC.m))
n_ref2_overlap <- sum(rownames(beta_pards) %in% rownames(centDHSbloodDMC.m))
cat("\nReference CpG overlap in the PARDS-only matrix:\n")
cat("  Ref1 (Epi/Fib/IC): ", n_ref1_overlap, "of", nrow(centEpiFibIC.m), "\n")
cat("  Ref2 (immune sub): ", n_ref2_overlap, "of", nrow(centDHSbloodDMC.m), "\n")

cat("\nRunning EpiDISH hepidish() with RPC on", n_pards, "PARDS specimens...\n")
hepidish_pards <- hepidish(
  beta.m   = beta_pards,
  ref1.m   = centEpiFibIC.m,
  ref2.m   = centDHSbloodDMC.m,
  h.CT.idx = 3,
  method   = "RPC"
)

cell_props_pards <- as.data.frame(hepidish_pards) %>%
  rownames_to_column("Sample_Name")

stopifnot(setequal(setdiff(names(cell_props_pards), "Sample_Name"), cell_cols))

cell_props_pards_out <- cell_props_pards %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID),
            by = "Sample_Name") %>%
  dplyr::select(Subject_ID, everything(), -Sample_Name)

cat("\nPARDS-only deconvolution estimates (rounded to 3dp):\n")
print(cell_props_pards_out %>% mutate(across(where(is.numeric), ~round(.x, 3))))
write_csv(cell_props_pards_out,
          file.path(SENS_TBL_DIR, "cell_composition_estimates_PARDS_only.csv"))

# -----------------------------------------------------------------------------
# S1.2 Deconvolution diagnostics and comparison with the full-cohort estimates
# -----------------------------------------------------------------------------
# Reviewer 3 comment 3.6 asks for the number of reference CpGs retained and the
# distributions and boundary values of the estimated fractions, so both are
# reported here for the PARDS-only fit.

prop_mat_pards <- as.matrix(
  cell_props_pards[match(pards_samples, cell_props_pards$Sample_Name),
                   cell_cols, drop = FALSE]
)
rownames(prop_mat_pards) <- pards_samples

prop_mat_full_sub <- as.matrix(
  cell_props[match(pards_samples, cell_props$Sample_Name),
             cell_cols, drop = FALSE]
)
rownames(prop_mat_full_sub) <- pards_samples

delta <- prop_mat_pards - prop_mat_full_sub
cat("\nDifference between PARDS-only and full-cohort fractions (same specimens):\n")
cat("  Max absolute difference:", signif(max(abs(delta)), 3), "\n")
cat("  Mean absolute difference:", signif(mean(abs(delta)), 3), "\n")
cat("  (EpiDISH deconvolves each specimen independently, so the estimates are\n",
    "   expected to be identical to numerical tolerance.)\n")

decon_summary <- tibble(
  cell_type  = cell_cols,
  mean_pct   = colMeans(prop_mat_pards) * 100,
  sd_pct     = apply(prop_mat_pards, 2, sd) * 100,
  median_pct = apply(prop_mat_pards, 2, median) * 100,
  min_pct    = apply(prop_mat_pards, 2, min) * 100,
  max_pct    = apply(prop_mat_pards, 2, max) * 100,
  n_at_zero  = colSums(prop_mat_pards <= 1e-8),
  n_at_one   = colSums(prop_mat_pards >= 1 - 1e-8),
  max_abs_diff_vs_full = apply(abs(delta), 2, max)
) %>% arrange(desc(mean_pct))

cat("\nPARDS-only fraction distributions and boundary counts:\n")
print(decon_summary %>% mutate(across(where(is.numeric), ~signif(.x, 3))))

row_sums <- rowSums(prop_mat_pards)
cat("\nRow sums of estimated fractions: min", signif(min(row_sums), 4),
    "max", signif(max(row_sums), 4), "\n")

save_support_table(
  decon_summary %>%
    transmute(`Cell type`  = cell_type,
              `Mean (%)`   = sprintf("%.1f", mean_pct),
              `SD (%)`     = sprintf("%.1f", sd_pct),
              `Median (%)` = sprintf("%.1f", median_pct),
              `Min (%)`    = sprintf("%.1f", min_pct),
              `Max (%)`    = sprintf("%.1f", max_pct),
              `n at 0`     = n_at_zero,
              `n at 1`     = n_at_one),
  file.path(SENS_TBL_DIR, "Supporting_deconvolution_PARDS_only"),
  caption = paste0("Supporting table. Cell-type fractions estimated by ",
                   "EpiDISH (hepidish, RPC) in the ", n_pards,
                   " PARDS specimens analyzed without control subjects."),
  note = paste0("Estimates use the two-stage centEpiFibIC and centDHSbloodDMC ",
                "references, with ", n_ref1_overlap, " and ", n_ref2_overlap,
                " reference CpGs respectively present in the filtered matrix. ",
                "Columns n at 0 and n at 1 count specimens whose estimated ",
                "fraction lies on the boundary of the simplex, where the ",
                "estimate is least reliable."))

# -----------------------------------------------------------------------------
# III.2 Deconvolution within the new-PARDS-only cohort
# -----------------------------------------------------------------------------

beta_new <- beta_collapsed[, new_pards, drop = FALSE]
m_new    <- m_collapsed[,   new_pards, drop = FALSE]

cat("\nRunning hepidish() on", n_new, "new PARDS specimens...\n")
hepidish_new <- hepidish(beta.m = beta_new,
                         ref1.m = centEpiFibIC.m,
                         ref2.m = centDHSbloodDMC.m,
                         h.CT.idx = 3, method = "RPC")

prop_mat_new <- as.matrix(hepidish_new[new_pards, cell_cols, drop = FALSE])
rownames(prop_mat_new) <- new_pards

# Estimates are per-specimen, so they must match the S1 values exactly.
delta_new <- prop_mat_new - prop_mat_pards[new_pards, cell_cols, drop = FALSE]
cat("Max |difference| vs the S1 PARDS-only estimates:",
    signif(max(abs(delta_new)), 3), "\n")

# -----------------------------------------------------------------------------
# III.3 Analytic sets, stage-one estimates, and display order
# -----------------------------------------------------------------------------
# ANALYTIC_SETS bundles the beta matrix, M-value matrix, and hepidish fractions
# for each cohort, and is the object every remaining section loops over.
#
# The stage-one (epithelial / fibroblast / immune) estimates are computed once
# here and reused twice: for the compositional log ratio in section III.7 and
# for the stage-one adjustment strategy in PART VI. hepidish partitions the
# stage-one immune fraction using stage two, so the immune total below equals
# the sum of the seven hepidish immune columns by construction; the identity is
# verified in PART VI.

prop_mat_full <- as.matrix(
  cell_props[match(colnames(beta_collapsed), cell_props$Sample_Name),
             cell_cols, drop = FALSE])
rownames(prop_mat_full) <- colnames(beta_collapsed)

ANALYTIC_SETS <- list(
  list(label = COHORT_LABELS[1], beta = beta_collapsed,
       m = m_collapsed,       prop = prop_mat_full),
  list(label = COHORT_LABELS[2], beta = beta_pards,
       m = m_collapsed_pards, prop = prop_mat_pards),
  list(label = COHORT_LABELS[3], beta = beta_new,
       m = m_new,             prop = prop_mat_new))

set_levels <- vapply(ANALYTIC_SETS, `[[`, character(1), "label")

# Cell types ordered by mean abundance in the full cohort. Used for every table
# that lists cell types, so the rows read the same way throughout.
ct_order <- names(sort(colMeans(prop_mat_full), decreasing = TRUE))

stage1_by_set <- lapply(ANALYTIC_SETS, function(s) {
  e1 <- epidish(beta.m = s$beta, ref.m = centEpiFibIC.m, method = "RPC")$estF
  e1[colnames(s$beta), , drop = FALSE]
})
names(stage1_by_set) <- set_levels

stage1_est_all <- map_dfr(set_levels, function(nm)
  as.data.frame(stage1_by_set[[nm]]) %>%
    rownames_to_column("Sample_Name") %>%
    mutate(analytic_set = nm))

cat("\nStage-one compartment means (%), by cohort:\n")
print(as.data.frame(stage1_est_all %>%
                      group_by(analytic_set) %>%
                      summarise(across(c(Epi, Fib, IC),
                                       ~round(mean(.x) * 100, 1)),
                                .groups = "drop")))


# -----------------------------------------------------------------------------
# III.4 Supplemental Table 4: reference CpGs retained
# -----------------------------------------------------------------------------

ref_overlap <- map_dfr(ANALYTIC_SETS, function(s) {
  tibble(analytic_set = s$label,
         n_specimens  = ncol(s$beta),
         probes_in_matrix = nrow(s$beta),
         ref1_total = nrow(centEpiFibIC.m),
         ref1_retained = sum(rownames(s$beta) %in% rownames(centEpiFibIC.m)),
         ref2_total = nrow(centDHSbloodDMC.m),
         ref2_retained = sum(rownames(s$beta) %in% rownames(centDHSbloodDMC.m)))
})
cat("\nReference CpG overlap:\n"); print(as.data.frame(ref_overlap))
write_csv(ref_overlap, file.path(SENS_TBL_DIR, "decon_reference_cpg_overlap.csv"))

# Column names are spelled out for the supplementary table because the raw
# names put ref1_retained (roughly 340) directly beside ref2_total (333), two
# numbers measuring different things, which reads as though both references
# were the same size.
save_pub_table(
  ref_overlap %>%
    transmute(
      `Analytic set`                   = analytic_set,
      `Specimens (n)`                  = n_specimens,
      `Probes in filtered matrix`      = probes_in_matrix,
      `Stage 1 reference CpGs (total)` = ref1_total,
      `Stage 1 retained (n)`           = ref1_retained,
      `Stage 1 retained (%)`           = sprintf("%.0f", 100 * ref1_retained / ref1_total),
      `Stage 2 reference CpGs (total)` = ref2_total,
      `Stage 2 retained (n)`           = ref2_retained,
      `Stage 2 retained (%)`           = sprintf("%.0f", 100 * ref2_retained / ref2_total)),
  file.path(TBL_DIR, "SuppTable4_reference_cpg_overlap"),
  caption = paste0("Supplemental Table 4. Reference CpGs retained on the ",
                   "EPIC v2 array for each stage of hierarchical ",
                   "deconvolution."),
  note = paste0("Stage 1 uses the centEpiFibIC reference to resolve ",
                "epithelial, fibroblast, and total immune compartments; ",
                "stage 2 uses the centDHSbloodDMC reference to resolve seven ",
                "immune subtypes within the stage 1 immune compartment. ",
                "Counts are identical across analytic sets because probe ",
                "filtering was performed at the cohort level. The stage 2 ",
                "reference is both smaller and less well retained, leaving ",
                "roughly 19 CpGs per immune subtype, which is the basis for ",
                "the stage 1 sensitivity analysis reported separately."))

# -----------------------------------------------------------------------------
# III.5 Table 2: fraction distributions and boundary values
# -----------------------------------------------------------------------------

frac_dist <- map_dfr(ANALYTIC_SETS, function(s) {
  p <- s$prop
  map_dfr(cell_cols, function(ct) {
    x <- p[, ct]
    tibble(analytic_set = s$label, cell_type = ct,
           mean_pct   = mean(x) * 100,
           median_pct = median(x) * 100,
           q25_pct    = quantile(x, 0.25) * 100,
           q75_pct    = quantile(x, 0.75) * 100,
           min_pct    = min(x) * 100,
           max_pct    = max(x) * 100,
           n_at_zero  = sum(x <= 1e-8),
           n_at_one   = sum(x >= 1 - 1e-8),
           pct_at_boundary = 100 * mean(x <= 1e-8 | x >= 1 - 1e-8))
  })
})
cat("\nFraction distributions and boundary counts:\n")
print(as.data.frame(frac_dist %>% mutate(across(where(is.numeric), ~signif(.x, 3)))))
write_csv(frac_dist, file.path(SENS_TBL_DIR, "decon_fraction_distributions.csv"))

rowsum_chk <- map_dfr(ANALYTIC_SETS, function(s) {
  rs <- rowSums(s$prop)
  tibble(analytic_set = s$label, min_rowsum = min(rs), max_rowsum = max(rs),
         max_abs_dev_from_1 = max(abs(rs - 1)))
})
cat("\nRow sums of estimated fractions:\n")
print(as.data.frame(rowsum_chk %>% mutate(across(where(is.numeric), ~signif(.x, 4)))))

table2_fractions <- frac_dist %>%
  mutate(analytic_set = factor(analytic_set, levels = set_levels),
         cell_type    = factor(cell_type,    levels = ct_order)) %>%
  arrange(analytic_set, cell_type) %>%
  transmute(Cohort            = as.character(analytic_set),
            `Cell type`       = as.character(cell_type),
            `Mean (%)`        = sprintf("%.1f", mean_pct),
            `Median (%)`      = sprintf("%.1f", median_pct),
            `IQR (%)`         = sprintf("%.1f \u2013 %.1f", q25_pct, q75_pct),
            `Range (%)`       = sprintf("%.1f \u2013 %.1f", min_pct, max_pct),
            `n at 0`          = as.character(n_at_zero),
            `n at 1`          = as.character(n_at_one),
            `At boundary (%)` = sprintf("%.1f", pct_at_boundary))

save_grouped_table(
  table2_fractions, "Cohort",
  file.path(TBL_DIR, "Table2_fraction_distributions"),
  caption = paste0("Table 2. Distributions and boundary values of estimated ",
                   "cell-type fractions."),
  note = paste0("Cell-type fractions were estimated by hierarchical EpiDISH ",
                "(hepidish) using robust partial correlation with the ",
                "centEpiFibIC and centDHSbloodDMC reference panels. Values are ",
                "percentages of the estimated composition. IQR is the ",
                "interquartile range and Range is the minimum to maximum ",
                "across specimens. Columns n at 0 and n at 1 count specimens ",
                "whose estimated fraction lies exactly on the boundary of the ",
                "simplex, where the estimate is least reliable; At boundary is ",
                "the percentage of specimens at either boundary. Cell types ",
                "are ordered by mean fraction in the full cohort. Row sums of ",
                "the estimated fractions deviated from one by no more than ",
                signif(max(rowsum_chk$max_abs_dev_from_1), 2), "."),
  widths = c(1.05, 0.75, 0.80, 1.00, 1.05, 0.55, 0.55, 0.90))


# -----------------------------------------------------------------------------
# III.6 Supplemental Tables 5 to 7: per-specimen fractions
# -----------------------------------------------------------------------------

# Reviewer 3 comment 3.6 asks for the distributions of all estimated fractions.
# The summary table (S6.2b) gives the distributions; these give the underlying
# per-specimen values. Three separate tables rather than one stacked table, so
# that no single supplementary table runs to 58 rows.
#
# Rows are ordered by epithelial fraction so the compositional gradient across
# specimens is visible by inspection; columns are ordered by mean abundance in
# the full cohort. A median row is appended to each table.

id_of <- setNames(sample_sheet$Subject_ID, sample_sheet$Sample_Name)

specimen_table <- function(prop_mat, cohort_label, file_stem, table_no) {
  d <- as.data.frame(prop_mat * 100)
  d <- d[, ct_order, drop = FALSE]           # mean-abundance order from S6.2b
  d$Subject <- unname(id_of[rownames(prop_mat)])
  d <- d[order(-d$Epi), , drop = FALSE]

  out <- tibble(Subject = d$Subject)
  for (ct in ct_order) out[[ct]] <- sprintf("%.1f", d[[ct]])

  med <- tibble(Subject = "Median")
  for (ct in ct_order) med[[ct]] <- sprintf("%.1f", median(d[[ct]]))
  out <- bind_rows(out, med)

  write_csv(out, file.path(TBL_DIR, paste0(file_stem, ".csv")))

  save_pub_table(
    out, file.path(TBL_DIR, file_stem),
    caption = paste0("Supplemental Table ", table_no,
                     ". Estimated cell-type fractions for each specimen: ",
                     cohort_label, "."),
    note = paste0("Values are percentages of the estimated cellular ",
                  "composition from hierarchical EpiDISH (hepidish) using ",
                  "robust partial correlation. Specimens are ordered by ",
                  "epithelial fraction and cell types by mean abundance in ",
                  "the full cohort. Fractions of 0.0 indicate an estimate at ",
                  "the boundary of the simplex rather than a rounded small ",
                  "value; boundary counts are given in the summary table. ",
                  "The final row gives the median across specimens in this ",
                  "cohort."))
  invisible(out)
}

cat("\n=== PART III.6. Per-specimen fraction tables ===\n")
tbl_full  <- specimen_table(prop_mat_full,  "full cohort (n = 24)",
                            "SuppTable5_fractions_by_specimen_full_cohort", "5")
tbl_pards <- specimen_table(prop_mat_pards, paste0("PARDS only (n = ", n_pards, ")"),
                            "SuppTable6_fractions_by_specimen_pards_only",  "6")
tbl_new   <- specimen_table(prop_mat_new,   paste0("new PARDS only (n = ", n_new, ")"),
                            "SuppTable7_fractions_by_specimen_new_pards",   "7")

cat("Wrote three per-specimen fraction tables to:", SENS_TBL_DIR, "\n")

# Range statement for the Results text
cat("\n--- For the Results text ---\n")
for (nm in names(list(`full cohort` = prop_mat_full,
                      `PARDS only`   = prop_mat_pards,
                      `new PARDS only` = prop_mat_new))) {
  p <- list(`full cohort` = prop_mat_full, `PARDS only` = prop_mat_pards,
            `new PARDS only` = prop_mat_new)[[nm]]
  cat(sprintf("%s: epithelial %.1f%% to %.1f%%, neutrophil %.1f%% to %.1f%%\n",
              nm, min(p[, "Epi"]) * 100, max(p[, "Epi"]) * 100,
              min(p[, "Neutro"]) * 100, max(p[, "Neutro"]) * 100))
}

# -----------------------------------------------------------------------------
# III.7 Supplemental Table 8: deconvolution goodness of fit
# -----------------------------------------------------------------------------

# hepidish is two-stage, so fit is assessed at each stage separately by
# reconstructing the observed beta values at the reference CpGs from the
# estimated fractions and the reference centroids:
#   stage 1: centEpiFibIC  x  (Epi, Fib, IC), IC = sum of the immune fractions
#   stage 2: centDHSbloodDMC x (immune fractions renormalized to sum to one)
# Per specimen we report the Pearson correlation and RMSE between observed and
# reconstructed beta. RPC is a robust regression, so a low correlation flags a
# specimen whose composition is poorly described by the reference panel.

immune_cols <- setdiff(cell_cols, c("Epi", "Fib"))

decon_fit <- map_dfr(ANALYTIC_SETS, function(s) {
  b <- s$beta; p <- s$prop

  # stage 1
  cp1  <- intersect(rownames(b), rownames(centEpiFibIC.m))
  R1   <- centEpiFibIC.m[cp1, , drop = FALSE]
  IC   <- rowSums(p[, immune_cols, drop = FALSE])
  F1   <- cbind(Epi = p[, "Epi"], Fib = p[, "Fib"], IC = IC)
  F1   <- F1[, colnames(R1), drop = FALSE]     # align to reference column order
  pred1 <- R1 %*% t(F1)
  obs1  <- b[cp1, rownames(F1), drop = FALSE]

  # stage 2
  cp2   <- intersect(rownames(b), rownames(centDHSbloodDMC.m))
  R2    <- centDHSbloodDMC.m[cp2, , drop = FALSE]
  F2raw <- p[, colnames(R2), drop = FALSE]
  denom <- rowSums(F2raw)
  F2    <- sweep(F2raw, 1, ifelse(denom > 0, denom, NA_real_), "/")
  pred2 <- R2 %*% t(F2)
  obs2  <- b[cp2, rownames(F2), drop = FALSE]

  tibble(
    analytic_set = s$label,
    Sample_Name  = colnames(b),
    stage1_r     = vapply(seq_len(ncol(b)), function(j)
                     suppressWarnings(cor(obs1[, j], pred1[, j])), numeric(1)),
    stage1_rmse  = vapply(seq_len(ncol(b)), function(j)
                     sqrt(mean((obs1[, j] - pred1[, j])^2)), numeric(1)),
    stage2_r     = vapply(seq_len(ncol(b)), function(j)
                     suppressWarnings(cor(obs2[, j], pred2[, j])), numeric(1)),
    stage2_rmse  = vapply(seq_len(ncol(b)), function(j)
                     sqrt(mean((obs2[, j] - pred2[, j])^2)), numeric(1)),
    immune_total = IC
  )
}) %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID), by = "Sample_Name") %>%
  dplyr::select(analytic_set, Subject_ID, everything(), -Sample_Name)

cat("\nDeconvolution goodness of fit, per specimen:\n")
print(as.data.frame(decon_fit %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
write_csv(decon_fit, file.path(SENS_TBL_DIR, "decon_goodness_of_fit.csv"))

fit_summary <- decon_fit %>%
  group_by(analytic_set) %>%
  summarise(across(c(stage1_r, stage1_rmse, stage2_r, stage2_rmse),
                   list(median = ~median(.x, na.rm = TRUE),
                        min    = ~min(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"),
            n_stage2_undefined = sum(is.na(stage2_r)), .groups = "drop")
cat("\nGoodness-of-fit summary by analytic set:\n")
print(as.data.frame(fit_summary %>% mutate(across(where(is.numeric), ~round(.x, 3)))))

save_pub_table(
  decon_fit %>%
    mutate(analytic_set = factor(analytic_set, levels = set_levels)) %>%
    arrange(analytic_set, desc(immune_total)) %>%
    transmute(`Analytic set`        = as.character(analytic_set),
              `Subject`             = Subject_ID,
              `Immune fraction (%)` = sprintf("%.1f", immune_total * 100),
              `Stage 1 correlation` = sprintf("%.3f", stage1_r),
              `Stage 2 correlation` = ifelse(is.na(stage2_r), "\u2014",
                                             sprintf("%.3f", stage2_r))),
  file.path(TBL_DIR, "SuppTable8_deconvolution_goodness_of_fit"),
  caption = paste0("Supplemental Table 8. Deconvolution goodness of fit for ",
                   "each specimen."),
  note = paste0("Expected beta values at the reference CpGs were computed as ",
                "the product of the reference centroids and each specimen's ",
                "estimated fractions, then correlated with the observed values ",
                "by Pearson correlation. Stage 1 uses the ",
                ref_overlap$ref1_retained[1], " retained centEpiFibIC CpGs and ",
                "stage 2 the ", ref_overlap$ref2_retained[1], " retained ",
                "centDHSbloodDMC CpGs. Stage 2 is undefined, shown as a dash, ",
                "for specimens whose immune compartment was estimated at zero. ",
                "Specimens are ordered within each analytic set by immune ",
                "fraction; stage 2 fit deteriorates as the immune compartment ",
                "shrinks."))


# -----------------------------------------------------------------------------
# III.8 Table 3: the compositional axis against clinical and technical variables
# -----------------------------------------------------------------------------

# The dominant axis of variation in these data is cell composition (S6.5). The
# reviewer's overadjustment concern assumes that axis reflects disease biology,
# specifically epithelial injury and neutrophil recruitment. An alternative is
# that it reflects variation in how the brushing was obtained: a single nasal
# brushing per subject, with no replicates, is subject to substantial sampling
# variation in depth, location, and mucus content.
#
# This block tests the compositional axis against a small prespecified set of
# clinical and technical variables. The immune-to-epithelial log ratio from the
# stage-one decomposition is used as the summary of composition, because the
# stage-one estimates are supported by the larger stage-one reference panel
# (see Supplemental Table 4) rather than the 131 retained stage-two CpGs.
#
# WHAT THIS CAN AND CANNOT SHOW
#   A clinical association would support the reviewer's biological reading.
#   A technical association would be direct evidence of a batch contribution.
#   No association with either is consistent with unstructured specimen-level
#   variation, i.e. sampling variability, but does not establish it: the
#   analysis is underpowered and the available clinical variables are coarse.
#   Report as exploratory and describe a null as "no correlate identified",
#   not as evidence of sampling bias.
#
# Five tests are prespecified per analytic set (three biological, two
# technical) and FDR corrected within set. Remaining variables are described
# without testing, to avoid adding the many small comparisons that Reviewers 1
# and 2 have objected to elsewhere.

cat("\n=== PART III.8. Compositional axis vs clinical and technical variables ===\n")

# -----------------------------------------------------------------------------
# S8.1 Resolve column names
# -----------------------------------------------------------------------------
pick1 <- function(cands) {
  hit <- cands[cands %in% names(sample_sheet)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

COLS_S8 <- list(
  severity = pick1(c("PARDS_Severity", "PARDS_Category", "Highest_PARDS_Category",
                     "PARDS_severity", "Severity", "OI_Category")),
  pelod    = pick1(c("PELOD2", "PELOD_2", "PELOD", "Highest_PELOD2", "PELOD2_max")),
  vfd      = pick1(c("VFD", "Ventilator_Free_Days", "VFD_28", "Vent_Free_Days")),
  slide    = pick1(c("Sentrix_ID", "Slide", "SentrixID", "Sentrix_Barcode")),
  position = pick1(c("Sentrix_Position", "Position", "Sentrix_Pos")),
  age      = pick1(c("Age", "Age_years", "age")),
  sex      = pick1(c("Sex", "sex", "Gender")),
  outcome  = pick1(c("Outcome", "Mortality", "Survival", "Vital_Status"))
)
cat("Resolved columns:\n"); print(unlist(COLS_S8))
if (any(is.na(unlist(COLS_S8))))
  cat("NOT FOUND (tests on these are skipped):",
      paste(names(COLS_S8)[is.na(unlist(COLS_S8))], collapse = ", "), "\n")

if (!exists("recode_severity")) {
  recode_severity <- function(x) {
    x <- trimws(as.character(x))
    out <- dplyr::case_when(
      is.na(x)                                   ~ NA_character_,
      grepl("^none|^no ", x, ignore.case = TRUE) ~ "None",
      grepl("^mild", x, ignore.case = TRUE)      ~ "Mild",
      grepl("^mod",  x, ignore.case = TRUE)      ~ "Moderate",
      grepl("^sev",  x, ignore.case = TRUE)      ~ "Severe",
      TRUE                                       ~ x)
    factor(out, levels = c("None", "Mild", "Moderate", "Severe"))
  }
}

# -----------------------------------------------------------------------------
# S8.2 Immune-to-epithelial log ratio from the stage-one decomposition
# -----------------------------------------------------------------------------
# The additive log ratio log(IC / Epi) is used rather than the raw ratio: the
# raw ratio is unbounded above and compresses everything below one, whereas the
# log ratio is symmetric about zero (equal immune and epithelial content).
# Because the fibroblast fraction is well under 1%, this is very close to the
# epithelial fraction on a logit scale.
#
# Two specimens have IC estimated at exactly zero, so a pseudocount is added.
# S8.5 repeats the analysis excluding them.

stage1_est <- stage1_est_all

PSEUDO <- 0.001   # 0.1 percentage point

alr_dat <- stage1_est %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID, Sample_Group,
                          all_of(unname(na.omit(unlist(COLS_S8))))),
            by = "Sample_Name") %>%
  mutate(
    analytic_set = factor(analytic_set, levels = set_levels),
    alr_immune_epi = log((IC + PSEUDO) / (Epi + PSEUDO)),
    ic_is_zero     = IC <= 1e-8,
    epi_is_one     = Epi >= 1 - 1e-8)

cat("\nImmune-to-epithelial log ratio, by analytic set:\n")
print(alr_dat %>% group_by(analytic_set) %>%
        summarise(n = n(),
                  median = median(alr_immune_epi),
                  min = min(alr_immune_epi), max = max(alr_immune_epi),
                  n_ic_zero = sum(ic_is_zero), .groups = "drop") %>%
        mutate(across(where(is.numeric), ~round(.x, 3))) %>% as.data.frame())

write_csv(alr_dat %>% dplyr::select(analytic_set, Subject_ID, Epi, Fib, IC,
                                    alr_immune_epi, ic_is_zero),
          file.path(SENS_TBL_DIR, "compositional_log_ratio.csv"))

# -----------------------------------------------------------------------------
# S8.3 Prespecified tests
# -----------------------------------------------------------------------------
# Continuous predictors: Spearman correlation, with a Fisher z approximate 95%
#   confidence interval.
# Categorical predictors: Kruskal-Wallis, with epsilon squared as the effect
#   size. Levels with fewer than 2 specimens are dropped before testing.

spearman_row <- function(y, x, label, domain, set_label) {
  ok <- !is.na(x) & !is.finite(x) == FALSE & !is.na(y)
  x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 5 || length(unique(x)) < 3)
    return(tibble(analytic_set = set_label, domain = domain, variable = label,
                  test = "Spearman", n = n, statistic = NA_real_,
                  effect = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                  p_value = NA_real_, note = "too few observations or levels"))
  ct  <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  se  <- 1 / sqrt(n - 3)
  tibble(analytic_set = set_label, domain = domain, variable = label,
         test = "Spearman", n = n, statistic = unname(ct$statistic),
         effect = rho,
         ci_low  = tanh(atanh(rho) - 1.96 * se),
         ci_high = tanh(atanh(rho) + 1.96 * se),
         p_value = ct$p.value, note = "")
}

kruskal_row <- function(y, g, label, domain, set_label, min_per_level = 2) {
  g <- as.character(g)
  keep <- !is.na(g) & !is.na(y)
  y <- y[keep]; g <- g[keep]
  tb <- table(g); good <- names(tb)[tb >= min_per_level]
  y <- y[g %in% good]; g <- g[g %in% good]
  k <- length(unique(g)); n <- length(y)
  if (k < 2 || n < 5)
    return(tibble(analytic_set = set_label, domain = domain, variable = label,
                  test = "Kruskal-Wallis", n = n, statistic = NA_real_,
                  effect = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                  p_value = NA_real_,
                  note = paste0("fewer than 2 testable levels (levels kept: ",
                                paste(good, collapse = ", "), ")")))
  kt  <- kruskal.test(y ~ factor(g))
  eps <- unname((kt$statistic - k + 1) / (n - k))
  tibble(analytic_set = set_label, domain = domain, variable = label,
         test = "Kruskal-Wallis", n = n, statistic = unname(kt$statistic),
         effect = eps, ci_low = NA_real_, ci_high = NA_real_,
         p_value = kt$p.value,
         note = paste0(k, " levels; sizes ",
                       paste(as.integer(tb[good]), collapse = "/")))
}

run_tests <- function(d, set_label) {
  y <- d$alr_immune_epi
  out <- list()

  if (!is.na(COLS_S8$severity))
    out[[length(out) + 1]] <- kruskal_row(y, recode_severity(d[[COLS_S8$severity]]),
                                          "PARDS severity", "Biological", set_label)
  if (!is.na(COLS_S8$pelod))
    out[[length(out) + 1]] <- spearman_row(y, suppressWarnings(as.numeric(d[[COLS_S8$pelod]])),
                                           "PELOD-2 score", "Biological", set_label)
  if (!is.na(COLS_S8$vfd))
    out[[length(out) + 1]] <- spearman_row(y, suppressWarnings(as.numeric(d[[COLS_S8$vfd]])),
                                           "Ventilator-free days", "Biological", set_label)
  if (!is.na(COLS_S8$slide))
    out[[length(out) + 1]] <- kruskal_row(y, d[[COLS_S8$slide]],
                                          "Array slide", "Technical", set_label)
  if (!is.na(COLS_S8$position))
    out[[length(out) + 1]] <- kruskal_row(y, d[[COLS_S8$position]],
                                          "Array position", "Technical", set_label)

  bind_rows(out) %>% mutate(p_fdr = p.adjust(p_value, method = "BH"))
}

s8_tests <- map_dfr(set_levels, function(sl)
  run_tests(alr_dat %>% filter(analytic_set == sl), sl))

cat("\nPrespecified tests of the compositional log ratio:\n")
print(as.data.frame(s8_tests %>% mutate(across(where(is.numeric), ~round(.x, 4)))))
write_csv(s8_tests, file.path(SENS_TBL_DIR, "compositional_axis_tests.csv"))

# Free comparison: PARDS vs control, full cohort only, descriptive
d_full <- alr_dat %>% filter(analytic_set == set_levels[1])
grp_desc <- d_full %>% group_by(Sample_Group) %>%
  summarise(n = n(), median_alr = median(alr_immune_epi),
            min = min(alr_immune_epi), max = max(alr_immune_epi), .groups = "drop")
cat("\nCompositional log ratio by group (full cohort, descriptive only):\n")
print(as.data.frame(grp_desc %>% mutate(across(where(is.numeric), ~round(.x, 3)))))

# -----------------------------------------------------------------------------
# S8.4 Descriptive summaries for the untested variables
# -----------------------------------------------------------------------------
desc_vars <- c(age = COLS_S8$age, sex = COLS_S8$sex, outcome = COLS_S8$outcome)
desc_vars <- desc_vars[!is.na(desc_vars)]

s8_desc <- map_dfr(names(desc_vars), function(nm) {
  col <- desc_vars[[nm]]
  d <- alr_dat %>% filter(analytic_set == set_levels[2])   # PARDS only
  v <- d[[col]]
  if (is.numeric(v) || !any(is.na(suppressWarnings(as.numeric(v))))) {
    vn <- suppressWarnings(as.numeric(v))
    tibble(variable = nm,
           summary = sprintf("Spearman rho = %.2f (descriptive, not tested)",
                             suppressWarnings(cor(vn, d$alr_immune_epi,
                                                  method = "spearman",
                                                  use = "complete.obs"))))
  } else {
    lv <- sort(unique(as.character(v[!is.na(v)])))
    tibble(variable = nm,
           summary = paste(vapply(lv, function(l)
             sprintf("%s: median %.2f (n = %d)", l,
                     median(d$alr_immune_epi[which(v == l)]), sum(v == l, na.rm = TRUE)),
             character(1)), collapse = "; "))
  }
})
cat("\nDescriptive only (PARDS-only set, no tests performed):\n")
print(as.data.frame(s8_desc))

# -----------------------------------------------------------------------------
# S8.5 Sensitivity: exclude the two zero-immune specimens
# -----------------------------------------------------------------------------
alr_nozero <- alr_dat %>% filter(!ic_is_zero)
s8_tests_nozero <- map_dfr(set_levels, function(sl)
  run_tests(alr_nozero %>% filter(analytic_set == sl), sl)) %>%
  mutate(analysis = "Excluding zero-immune specimens")

cat("\nSensitivity, excluding specimens with immune fraction estimated at zero:\n")
print(as.data.frame(s8_tests_nozero %>% mutate(across(where(is.numeric), ~round(.x, 4)))))
write_csv(s8_tests_nozero,
          file.path(SENS_TBL_DIR, "compositional_axis_tests_excluding_zeros.csv"))

# -----------------------------------------------------------------------------
# S8.6 Publication table
# -----------------------------------------------------------------------------
save_pub_table(
  s8_tests %>%
    transmute(`Analytic set` = as.character(analytic_set),
              `Domain`       = domain,
              `Variable`     = variable,
              `Test`         = test,
              `n`            = n,
              `Effect size`  = ifelse(is.na(effect), "\u2014",
                                      sprintf("%.2f", effect)),
              `95% CI`       = ifelse(is.na(ci_low), "\u2014",
                                      sprintf("%.2f to %.2f", ci_low, ci_high)),
              `p`            = ifelse(is.na(p_value), "\u2014",
                                      sprintf("%.3f", p_value)),
              `FDR p`        = ifelse(is.na(p_fdr), "\u2014",
                                      sprintf("%.3f", p_fdr))),
  file.path(TBL_DIR, "Table3_compositional_axis_tests"),
  caption = paste0("Table 3. Association of the immune-to-",
                   "epithelial compositional log ratio with prespecified ",
                   "clinical and technical variables."),
  note = paste0("The compositional summary is the additive log ratio ",
                "log(immune / epithelial) from the stage-one deconvolution, ",
                "with a pseudocount of ", PSEUDO, " added to both fractions ",
                "because two specimens had the immune fraction estimated at ",
                "zero. Effect sizes are Spearman rho for continuous variables ",
                "(with an approximate Fisher z confidence interval) and ",
                "epsilon squared for Kruskal-Wallis tests. Five tests were ",
                "prespecified per analytic set, three biological and two ",
                "technical, with Benjamini-Hochberg correction applied within ",
                "each set. Age, sex, comorbidity, acute condition, and outcome ",
                "were summarized descriptively without testing. These analyses ",
                "are exploratory and are not powered to exclude moderate ",
                "associations."))

# -----------------------------------------------------------------------------
# S8.7 Figures
# -----------------------------------------------------------------------------
# No figure is written for this section; Table 3 reports every test.

# -----------------------------------------------------------------------------
# S8.8 Checkpoint and text-ready summary
# -----------------------------------------------------------------------------
saveRDS(list(alr_dat = alr_dat, tests = s8_tests,
             tests_nozero = s8_tests_nozero, group_desc = grp_desc,
             descriptive = s8_desc, cols = COLS_S8, pseudocount = PSEUDO),
        file.path(RDS_DIR, "checkpoint_compositional_axis_S8.rds"))

cat("\n--- Summary for the response letter ---\n")
sig <- s8_tests %>% filter(!is.na(p_fdr), p_fdr < 0.05)
cat(sprintf("Prespecified tests run: %d across %d analytic sets.\n",
            nrow(s8_tests), length(set_levels)))
if (nrow(sig) == 0) {
  cat("No variable was associated with the compositional log ratio after FDR\n")
  cat("correction, in any analytic set. Report as: no clinical or technical\n")
  cat("correlate of the compositional axis was identified. Do NOT state that\n")
  cat("this demonstrates sampling variation; it is consistent with it.\n")
} else {
  cat("Associations surviving FDR correction:\n")
  print(as.data.frame(sig %>% dplyr::select(analytic_set, domain, variable,
                                            effect, p_value, p_fdr)))
  cat("\nIf the surviving association is technical, this is direct evidence of\n")
  cat("a batch contribution. If biological, it supports the reviewer's reading\n")
  cat("that the compositional axis reflects disease processes.\n")
}

cat("\nPART III.8 complete:", format(Sys.time()), "\n")


# =============================================================================

# -----------------------------------------------------------------------------
# III.9 Table 4: unadjusted PC1 against each estimated cell fraction
# -----------------------------------------------------------------------------

# Establishes that the dominant axis of the unadjusted methylome is cellular
# composition. PC1 is taken from the unadjusted PCA in each cohort and
# correlated with every estimated cell fraction.
#
# SIGN CONVENTION: the direction of a principal component is arbitrary and PC1
# emerged with opposite orientation in the 13-specimen cohort. PC1 is therefore
# oriented in each cohort so that it correlates positively with the epithelial
# fraction, which makes the rows comparable across cohorts and the neutrophil
# correlation consistently negative. This is also the value to quote for the
# neutrophil-PC1 sign raised in review.

cat("\n=== PART III.9. PC1 vs cell composition (Table 4) ===\n")

pc1_composition <- map_dfr(ANALYTIC_SETS, function(s) {
  pv    <- rowVars(s$m, na.rm = TRUE)
  m_top <- s$m[pv >= quantile(pv, 0.95, na.rm = TRUE), , drop = FALSE]
  pca   <- prcomp(t(m_top), center = TRUE, scale. = FALSE)
  ve    <- (pca$sdev^2) / sum(pca$sdev^2)

  pc1  <- pca$x[, 1]
  epi  <- s$prop[names(pc1), "Epi"]
  flip <- sign(cor(pc1, epi))
  if (is.na(flip) || flip == 0) flip <- 1
  pc1 <- pc1 * flip

  map_dfr(cell_cols, function(ct) {
    x <- s$prop[names(pc1), ct]
    if (sd(x) == 0)
      return(tibble(cohort = s$label, n = ncol(s$m), pc1_pct = ve[1] * 100,
                    cell_type = ct, r = NA_real_, ci_low = NA_real_,
                    ci_high = NA_real_, p_value = NA_real_,
                    note = "no variation in this cohort"))
    ctest <- cor.test(pc1, x)
    tibble(cohort = s$label, n = ncol(s$m), pc1_pct = ve[1] * 100,
           cell_type = ct, r = unname(ctest$estimate),
           ci_low = ctest$conf.int[1], ci_high = ctest$conf.int[2],
           p_value = ctest$p.value, note = "")
  })
}) %>%
  group_by(cohort) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(cohort    = factor(cohort, levels = set_levels),
         cell_type = factor(cell_type, levels = ct_order)) %>%
  arrange(cohort, cell_type)

cat("\nPearson correlation of unadjusted PC1 with each cell fraction:\n")
print(as.data.frame(pc1_composition %>%
                      mutate(across(where(is.numeric), ~round(.x, 4)))))
write_csv(pc1_composition, file.path(SENS_TBL_DIR, "pc1_vs_composition.csv"))

save_pub_table(
  pc1_composition %>%
    transmute(`Cohort`     = as.character(cohort),
              `Cell type`  = as.character(cell_type),
              `r`          = ifelse(is.na(r), "\u2014", sprintf("%.2f", r)),
              `95% CI`     = ifelse(is.na(ci_low), "\u2014",
                                    sprintf("%.2f to %.2f", ci_low, ci_high)),
              `p`          = ifelse(is.na(p_value), "\u2014",
                                    format.pval(p_value, digits = 2, eps = 1e-4)),
              `FDR p`      = ifelse(is.na(p_fdr), "\u2014",
                                    format.pval(p_fdr, digits = 2, eps = 1e-4))),
  file.path(TBL_DIR, "Table4_pc1_vs_composition"),
  caption = paste0("Table 4. Correlation of the first principal ",
                   "component of the unadjusted methylome with each estimated ",
                   "cell-type fraction."),
  note = paste0("Values are Pearson correlations between PC1 of the ",
                "unadjusted methylome and the estimated fraction of each cell ",
                "type, computed within each cohort. Because the direction of a ",
                "principal component is arbitrary, PC1 was oriented in each ",
                "cohort so that it correlates positively with the epithelial ",
                "fraction; the neutrophil correlation is negative throughout ",
                "because these two fractions vary inversely. P values are ",
                "Benjamini-Hochberg corrected within cohort across the ",
                length(cell_cols), " cell types tested. Cell types are ordered ",
                "by mean abundance in the full cohort."))

cat("\n--- For the Results text ---\n")
for (sl in set_levels) {
  d <- pc1_composition %>% filter(cohort == sl)
  e <- d %>% filter(cell_type == "Epi")
  nn <- d %>% filter(cell_type == "Neutro")
  cat(sprintf("%s (PC1 = %.1f%%): epithelial r = %.2f (p %s), neutrophil r = %.2f (p %s)\n",
              sl, d$pc1_pct[1], e$r, format.pval(e$p_value, digits = 2, eps = 1e-4),
              nn$r, format.pval(nn$p_value, digits = 2, eps = 1e-4)))
}
sig <- pc1_composition %>% dplyr::filter(!is.na(p_fdr), p_fdr < 0.05) %>%
  dplyr::count(cell_type, name = "cohorts_significant")
cat("\nCell types correlated with PC1 after FDR correction, by number of cohorts:\n")
print(as.data.frame(sig))

saveRDS(list(ref_overlap = ref_overlap, frac_dist = frac_dist,
             rowsum_chk = rowsum_chk, decon_fit = decon_fit,
             fit_summary = fit_summary, pc1_composition = pc1_composition,
             stage1_by_set = stage1_by_set, alr_tests = s8_tests),
        file.path(RDS_DIR, "checkpoint_deconvolution_PART_III.rds"))

cat("\nPART III complete:", format(Sys.time()), "\n")


# =============================================================================


# =============================================================================
# PART IV. COMPOSITION ADJUSTMENT, PCA, AND CLUSTERING
# =============================================================================


# -----------------------------------------------------------------------------
# IV.1 Full cohort: PCA of the adjusted methylome (Figure 2)
# -----------------------------------------------------------------------------

# from the raw analysis because the composition-driven probes won't dominate).

probe_var_adj <- rowVars(m_adjusted, na.rm = TRUE)
top5_adj_idx  <- which(probe_var_adj >= quantile(probe_var_adj, 0.95, na.rm = TRUE))
m_adj_top5    <- m_adjusted[top5_adj_idx, ]
cat("\nTop 5% variable CpGs (adjusted):", length(top5_adj_idx), "\n")

pca_adj <- prcomp(t(m_adj_top5), center = TRUE, scale. = FALSE)
var_exp_adj <- (pca_adj$sdev^2) / sum(pca_adj$sdev^2)

variance_compare <- tibble(
  PC          = paste0("PC", 1:10),
  raw_pct     = round(var_exp_raw[1:10] * 100, 1),
  adjusted_pct = round(var_exp_adj[1:10] * 100, 1)
)
cat("\nVariance explained by first 10 PCs (raw vs adjusted):\n")
print(variance_compare)
write_csv(variance_compare, file.path(TBL_DIR, "pca_variance_comparison.csv"))


pca_adj_df <- as.data.frame(pca_adj$x) %>%
  rownames_to_column("Sample_Name") %>%
  left_join(sample_sheet, by = "Sample_Name") %>%
  mutate(label = fig_label(Subject_ID))

# Targeted nudge: pull the ARDS026 label up and to the left, clear of its own
# dot in the upper-left corner. min.segment.length = 0 keeps a leader stem from
# the dot to the moved label. Tune the magnitudes (or flip the signs) to taste.
# The panel expansion gives the moved label room so it is not clipped at the edge.
# Per-label manual nudges for THIS layout. ggrepel won't reroute a stem that
# crosses a point, so we move offending labels by hand. Keyed by label; any
# label not listed gets 0. Tune magnitudes / flip signs to taste.
# Keyed by the rendered label, so they follow STRIP_UNDERSCORE.
nudge_x_map <- setNames(c(-12,  16, -10,  10),
                        fig_label(c("ARDS_026", "ARDS_004",
                                    "Control_001", "Control_022")))
nudge_y_map <- setNames(c( 12, -10,   8,  -8),
                        fig_label(c("ARDS_026", "ARDS_004",
                                    "Control_001", "Control_022")))

lab_nudge_x <- unname(nudge_x_map[pca_adj_df$label]); lab_nudge_x[is.na(lab_nudge_x)] <- 0
lab_nudge_y <- unname(nudge_y_map[pca_adj_df$label]); lab_nudge_y[is.na(lab_nudge_y)] <- 0

p_pca <- ggplot(pca_adj_df, aes(x = PC1, y = PC2,
                                color = Sample_Group, label = label)) +
  geom_point(size = 3) +
  geom_text_repel(size = 2.7, show.legend = FALSE,
                  min.segment.length = 0,
                  segment.size = 0.3, segment.color = "grey55",
                  box.padding = 0.5, point.padding = 0.25,
                  nudge_x = lab_nudge_x, nudge_y = lab_nudge_y,
                  max.overlaps = Inf, seed = 1) +
  scale_color_manual(values = c(Control = FIG1_BLUE, ARDS = FIG1_ORANGE),
                     name = "Group") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(x = paste0("PC1 (", round(var_exp_adj[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_exp_adj[2] * 100, 1), "%)"),
       title = paste0("Figure 2. PCA of Adjusted Methylome, Full Cohort (n = ",
                      nrow(pca_adj_df), ")")) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(FIG_DIR, "Figure2_PCA_adjusted.pdf"),
       p_pca, width = 7, height = 5, device = cairo_pdf, dpi = 600)
ggsave(file.path(FIG_DIR, "Figure2_PCA_adjusted.png"),
       p_pca, width = 7, height = 5, dpi = 600)
cat("Figure 2 written: Figure2_PCA_adjusted.pdf / .png\n")

# -----------------------------------------------------------------------------
# IV.2 Full cohort: hierarchical clustering and silhouette (Figure 3)
# -----------------------------------------------------------------------------

# Don't assume k=2. Test k = 2..6 and pick by mean silhouette width.
# Silhouette > 0.5 = strong structure; 0.25-0.5 = weak; < 0.25 = essentially
# no structure (samples form a continuum).

sample_dist_adj <- dist(t(m_adj_top5), method = "euclidean")
sample_hc_adj   <- hclust(sample_dist_adj, method = "ward.D2")

sil_table <- map_dfr(2:6, function(k) {
  cl  <- cutree(sample_hc_adj, k = k)
  if (length(unique(cl)) < k) return(tibble(k = k, mean_silhouette = NA_real_))
  sil <- silhouette(cl, sample_dist_adj)
  tibble(k = k, mean_silhouette = mean(sil[, 3]))
})

cat("\nSilhouette width by k (composition-adjusted clustering):\n")
print(sil_table %>% mutate(mean_silhouette = round(mean_silhouette, 3)))
write_csv(sil_table, file.path(TBL_DIR, "silhouette_by_k.csv"))

best_k <- sil_table$k[which.max(sil_table$mean_silhouette)]
best_sil <- max(sil_table$mean_silhouette, na.rm = TRUE)
cat("Best k by silhouette:", best_k, "(width =", round(best_sil, 3), ")\n")
cat("Interpretation: > 0.5 strong, 0.25-0.5 weak, < 0.25 essentially absent\n")

# Use best_k, with fallback to k=2 for downstream interpretability
chosen_k <- if (best_sil >= 0.25) best_k else 2
cat("Using k =", chosen_k, "for downstream analyses\n")

group_adj <- paste0("AdjGroup", cutree(sample_hc_adj, k = chosen_k))
sample_sheet$Adjusted_Group <- group_adj[match(sample_sheet$Sample_Name,
                                               names(cutree(sample_hc_adj, k = chosen_k)))]

cat("\nAdjusted cluster assignments:\n")
print(table(sample_sheet$Adjusted_Group,
            sample_sheet$Sample_Group, dnn = c("AdjustedGroup", "Sample_Group")))

write_csv(sample_sheet %>% dplyr::select(Subject_ID, Sample_Name, Sample_Group,
                                         Adjusted_Group, everything()),
          file.path(TBL_DIR, "sample_sheet_with_adjusted_groups.csv"))
# -----------------------------------------------------------------------------
# Annotated dendrogram: cluster assignments (bottom strip)
# -----------------------------------------------------------------------------
# Honest depiction of the (weak) k = chosen_k partition. The tree is left BLACK
# so branch color does not overstate the separation; the k = chosen_k cut is
# drawn; and a single strip beneath the leaves encodes the cluster assignment
# stored in Adjusted_Group. Diagnosis is readable from the Subject_ID labels.
# Pairs with Table 5 (width 0.233 at k = 2, below the 0.25
# floor): the figure shows the split, the table shows it is shallow.
if (!requireNamespace("dendextend", quietly = TRUE))
  stop("Install dendextend for this figure: install.packages('dendextend')")
suppressPackageStartupMessages(library(dendextend))

# (a) Dendrogram; leaves are Sample_Names at this stage
dend <- as.dendrogram(sample_hc_adj)
leaf_samples <- labels(dend)          # leaf (plotting) order, Sample_Name terms

# (b) Cluster assignment in leaf order (named by Sample_Name)
adj_assign    <- cutree(sample_hc_adj, k = chosen_k)
clusters_leaf <- adj_assign[leaf_samples]

# (c) Cluster labels and colours, shared with Figures 4 and 5
lab_full        <- cluster_labels(adj_assign)
full_group_cols <- setNames(CLUSTER_PALETTE[seq_along(lab_full)],
                            unname(lab_full))

# (d) Relabel leaves to Subject_ID for display (leaf order unchanged)
labels(dend) <- fig_label(unname(setNames(sample_sheet$Subject_ID,
                                          sample_sheet$Sample_Name)[labels(dend)]))

# (e) Cut height for the k-cluster solution (midway between bracketing merges)
n_leaves <- length(leaf_samples)
hts      <- sort(sample_hc_adj$height)
cut_h    <- mean(c(hts[n_leaves - chosen_k], hts[n_leaves - chosen_k + 1]))

# (f) Rebuilt as a ggplot so Figures 3, 4B and 5B share one visual grammar.
#     The base-R version used a different typeface, rotated y-axis tick labels
#     and colored_bars(), which made Figure 3 read as though it came from a
#     different paper. The tree is still left BLACK so branch colour does not
#     overstate the separation; the strip beneath the leaves carries the
#     assignment. Pairs with Table 5 (0.233 at k = 2, below the 0.25 floor):
#     the figure shows the split, the table shows it is shallow.
ggd_f  <- dendextend::as.ggdend(dend)
segs_f <- ggd_f$segments
labs_f <- ggd_f$labels %>%
  mutate(group_lab = unname(lab_full[as.character(clusters_leaf)]))

max_h_f  <- max(segs_f$y, na.rm = TRUE)
strip_hf <- max_h_f * 0.035
strip_yf <- -strip_hf / 2 - max_h_f * 0.01
label_yf <- -strip_hf * 1.5

fig3_dendro <- ggplot() +
  geom_segment(data = segs_f, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.35, color = "black") +
  geom_hline(yintercept = cut_h, linetype = "dashed", color = "grey40") +
  geom_tile(data = labs_f, aes(x = x, y = strip_yf, fill = group_lab),
            height = strip_hf, width = 0.85) +
  geom_text(data = labs_f, aes(x = x, y = label_yf, label = label),
            angle = 90, hjust = 1, vjust = 0.5, size = 2.4) +
  scale_fill_manual(values = full_group_cols,
                    name = paste0("Cluster (k = ", chosen_k, ")")) +
  scale_x_continuous(expand = expansion(add = 0.7)) +
  scale_y_continuous(breaks = pretty(c(0, max_h_f))) +
  coord_cartesian(ylim = c(-max_h_f * 0.45, max_h_f * 1.02), clip = "off") +
  # theme_classic would draw the y axis line down past 0 alongside the leaf
  # labels; draw it by hand so it stops at the top of the tree, as the base-R
  # version did.
  annotate("segment", x = -Inf, xend = -Inf, y = 0, yend = max_h_f,
           linewidth = 0.5) +
  labs(x = NULL, y = "Height",
       title = paste0("Figure 3. Hierarchical Clustering of the Adjusted ",
                      "Methylome, Full Cohort (n = ", n_leaves,
                      ", k = ", chosen_k, ")")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x         = element_blank(),
        axis.ticks.x        = element_blank(),
        axis.line.x         = element_blank(),
        axis.line.y         = element_blank(),
        legend.position     = "bottom",
        plot.title          = element_text(hjust = 0.5, face = "bold"),
        plot.title.position = "plot")

ggsave(file.path(FIG_DIR, "Figure3_dendrogram_adjusted.pdf"), fig3_dendro,
       width = 11, height = 6, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure3_dendrogram_adjusted.png"), fig3_dendro,
       width = 11, height = 6, dpi = 600)
cat("Figure 3 written: Figure3_dendrogram_adjusted.pdf / .png\n")

# -----------------------------------------------------------------------------
# IV.3 PARDS only: adjustment re-derived within the cohort
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# The two reductions applied to the covariate set are exactly those stated in
# the manuscript methods, applied unchanged to the PARDS-only estimates:
#   (i)  the cell type with the largest mean proportion is dropped to serve as
#        the reference and break the sum-to-one collinearity;
#   (ii) any remaining cell type proportion with a variance below 1e-4 across
#        samples is dropped, these being present at only low proportions and
#        not contributing informative adjustment.
# Both are recomputed inside the subset rather than inherited, because the
# largest-mean column and the set of low-variance columns can differ once the
# controls are removed. Variances are computed on the proportion scale (0 to 1),
# matching the scale on which the 1e-4 threshold was defined.

VAR_THRESHOLD <- 1e-4

drop_col_pards <- cell_cols[which.max(colMeans(prop_mat_pards))]
cat("\nDropping reference cell type (largest mean within PARDS):",
    drop_col_pards, "(", round(colMeans(prop_mat_pards)[drop_col_pards], 3), ")\n")
if (exists("prop_for_adj")) {
  cat("  Full-cohort adjustment retained:",
      paste(colnames(prop_for_adj), collapse = ", "), "\n")
}

candidate_pards <- setdiff(cell_cols, drop_col_pards)
col_vars_pards  <- apply(prop_mat_pards[, candidate_pards, drop = FALSE], 2, var)
nzv_pards       <- names(col_vars_pards)[col_vars_pards < VAR_THRESHOLD]
if (length(nzv_pards) > 0)
  cat("  Dropping low-variance composition columns:",
      paste(nzv_pards, collapse = ", "), "\n")

keep_pards     <- setdiff(candidate_pards, nzv_pards)
prop_adj_pards <- prop_mat_pards[, keep_pards, drop = FALSE]
cat("  Composition covariates retained:", paste(keep_pards, collapse = ", "),
    "(", ncol(prop_adj_pards), "columns )\n")
cat("  Residual degrees of freedom:", n_pards - ncol(prop_adj_pards) - 1, "\n")

# Covariate disposition table. Reporting only: the selection rules above are
# unchanged. Records the mean and variance behind every keep/drop decision so
# the covariate set used during residualization is fully specified, including
# any fraction whose variance falls close to the threshold.
covariate_table <- tibble(
  cell_type = cell_cols,
  mean_prop = colMeans(prop_mat_pards)[cell_cols],
  variance  = apply(prop_mat_pards[, cell_cols, drop = FALSE], 2, var),
  n_nonzero = colSums(prop_mat_pards[, cell_cols, drop = FALSE] > 1e-8)
) %>%
  mutate(
    disposition = case_when(
      cell_type == drop_col_pards ~ "Dropped: reference (largest mean)",
      cell_type %in% nzv_pards    ~ paste0("Dropped: variance < ", format(VAR_THRESHOLD)),
      TRUE                        ~ "Retained as covariate"),
    ratio_to_threshold = variance / VAR_THRESHOLD
  ) %>%
  arrange(desc(mean_prop))

cat("\nCovariate disposition (composition adjustment, PARDS only):\n")
print(covariate_table %>%
        mutate(mean_pct = signif(mean_prop * 100, 3),
               variance = signif(variance, 3),
               ratio_to_threshold = signif(ratio_to_threshold, 3)) %>%
        dplyr::select(cell_type, mean_pct, variance, ratio_to_threshold,
                      n_nonzero, disposition) %>%
        as.data.frame())

write_csv(covariate_table,
          file.path(SENS_TBL_DIR, "covariate_disposition_PARDS_only.csv"))

save_support_table(
  covariate_table %>%
    transmute(`Cell type`             = cell_type,
              `Mean proportion (%)`   = sprintf("%.2f", mean_prop * 100),
              `Variance`              = formatC(variance, format = "e", digits = 2),
              `Specimens nonzero (n)` = n_nonzero,
              `Disposition`           = disposition),
  file.path(SENS_TBL_DIR, "Supporting_covariates_PARDS_only"),
  caption = paste0("Supporting table. Cell-type proportions and their ",
                   "disposition in the composition adjustment, PARDS ",
                   "participants only (n = ", n_pards, ")."),
  note = paste0("Cell-type proportions were supplied to limma::removeBatchEffect ",
                "as continuous covariates. Two reductions were applied to avoid ",
                "rank deficiency: the cell type with the largest mean proportion ",
                "was dropped to serve as the reference and break the sum-to-one ",
                "collinearity, and any remaining proportion with a variance below ",
                "1e-4 across specimens was dropped. Variances are computed on the ",
                "proportion scale. The adjustment retained ",
                ncol(prop_adj_pards), " covariates, leaving ",
                n_pards - ncol(prop_adj_pards) - 1, " residual degrees of freedom."))

m_adj_pards <- removeBatchEffect(m_collapsed_pards, covariates = prop_adj_pards)

ok_pards <- rowSums(is.na(m_adj_pards) | !is.finite(m_adj_pards)) == 0
if (sum(!ok_pards) > 0)
  cat("  Dropping", sum(!ok_pards), "probes with NAs introduced by removeBatchEffect\n")
m_adj_pards <- m_adj_pards[ok_pards, , drop = FALSE]
cat("  PARDS-only analytic matrix:", nrow(m_adj_pards), "probes x",
    ncol(m_adj_pards), "samples\n")

saveRDS(m_adj_pards, file.path(RDS_DIR, "m_adjusted_composition_PARDS_only.rds"))

# -----------------------------------------------------------------------------
# IV.4 PARDS only: PCA, clustering, silhouette
# -----------------------------------------------------------------------------

res_pards <- cluster_diagnostics(m_adj_pards, "PARDS only")

# sil_table is the full-cohort silhouette computed in section IV.2; it is used
# directly rather than re-read from disk so the two can never disagree.
sil_full <- sil_table %>%
  mutate(analysis = paste0("Full cohort (n = ", ncol(m_adjusted),
                           ", controls included)"))

sil_compare <- bind_rows(
  sil_full,
  res_pards$sil %>% mutate(analysis = paste0("PARDS only (n = ", n_pards, ")"))
) %>%
  dplyr::select(analysis, k, mean_silhouette)

write_csv(sil_compare, file.path(SENS_TBL_DIR, "silhouette_by_k_PARDS_only.csv"))

suppTable_sil_sens <- sil_compare %>%
  mutate(mean_silhouette = sprintf("%.3f", mean_silhouette)) %>%
  pivot_wider(names_from = analysis, values_from = mean_silhouette) %>%
  rename(`Number of clusters (k)` = k)

save_support_table(
  suppTable_sil_sens,
  file.path(SENS_TBL_DIR, "Supporting_silhouette_PARDS_only"),
  caption = paste0("Supporting table. Mean silhouette width for k = 2 to 6 ",
                   "in the full cohort and in the PARDS-only sensitivity ",
                   "analysis."),
  note = paste0("Cell-type deconvolution, composition adjustment, and the top ",
                "5% variance filter were all re-derived using the ", n_pards,
                " PARDS specimens alone. Silhouette widths were computed on ",
                "Ward.D2 hierarchical clustering of euclidean distances. ",
                "Interpretation thresholds: greater than 0.5 strong, 0.25 to ",
                "0.5 weak, less than 0.25 essentially absent."))

# -----------------------------------------------------------------------------
# IV.5 PARDS only: concordance with the full-cohort partition
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Does dropping the controls change who clusters with whom? Compare the
# PARDS-only cut at the same k used in the primary analysis against the
# Adjusted_Group assignment from the 24-sample run, restricted to PARDS.

k_compare <- chosen_k

# Cut once and reuse. cutree() returns a vector named by Sample_Name, ordered
# to match the rows of res_pards$dist, which is what silhouette() indexes below.
cl_pards <- cutree(res_pards$hc, k = k_compare)
stopifnot(!is.null(names(cl_pards)),
          identical(names(cl_pards), labels(res_pards$dist)))

sil_at_k <- res_pards$sil$mean_silhouette[res_pards$sil$k == k_compare]

assign_pards <- tibble(
  Sample_Name      = names(cl_pards),
  PARDS_only_Group = paste0("PardsGroup", unname(cl_pards))
)

concordance <- sample_sheet %>%
  filter(Sample_Name %in% pards_samples) %>%
  dplyr::select(Sample_Name, Subject_ID, Adjusted_Group) %>%
  left_join(assign_pards, by = "Sample_Name") %>%
  arrange(Subject_ID)

cat("\nCluster assignment at k =", k_compare, ": full-cohort vs PARDS-only\n")
print(table(concordance$Adjusted_Group, concordance$PARDS_only_Group,
            dnn = c("Full cohort", "PARDS only")))

ari <- adjusted_rand_index(concordance$Adjusted_Group,
                           concordance$PARDS_only_Group)
cat("Adjusted Rand index, full cohort vs PARDS-only:", round(ari, 3), "\n")

write_csv(concordance,
          file.path(SENS_TBL_DIR, "cluster_assignment_concordance_PARDS_only.csv"))

# Per-sample silhouette at k_compare. Negative values flag specimens that sit
# closer to the other cluster than to their own. silhouette() returns a matrix
# with no rownames, but its rows are in the order of the clustering vector
# passed in, so the sample names come from names(cl_pards).
sil_obj <- silhouette(cl_pards, res_pards$dist)
stopifnot(nrow(sil_obj) == length(cl_pards))
sil_sample <- tibble(
  Sample_Name = names(cl_pards),
  cluster     = paste0("PardsGroup", sil_obj[, 1]),
  neighbor    = paste0("PardsGroup", sil_obj[, 2]),
  sil_width   = round(sil_obj[, 3], 3)
) %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID), by = "Sample_Name") %>%
  dplyr::select(Subject_ID, cluster, neighbor, sil_width) %>%
  arrange(cluster, desc(sil_width))

cat("\nPer-sample silhouette widths (k =", k_compare, "):\n")
print(as.data.frame(sil_sample))
write_csv(sil_sample,
          file.path(SENS_TBL_DIR, "per_sample_silhouette_PARDS_only.csv"))


# -----------------------------------------------------------------------------
# IV.6 PARDS only: figures (Figure 4)
# -----------------------------------------------------------------------------

# Publication settings: cairo_pdf for vector, 600 dpi PNG for submission.
# Only figures cited in the manuscript are written.

# (No silhouette line plot is written; Table 5 reports silhouette width and
#  cluster size for k = 2 to 6 in all three cohorts.)

# -----------------------------------------------------------------------------
# S1.7 Paneled figure: PCA (A) + dendrogram (B)
# -----------------------------------------------------------------------------
# Panel B is rebuilt as a ggplot via dendextend::as.ggdend so both panels are
# grid objects and patchwork can align them. Both panels map the cluster to
# fill and share one scale, but patchwork only MERGES guides whose key glyphs
# match, and panel A's points draw a point key while panel B's tiles draw a
# rect key -- which is why two identical legends used to print. Panel B's guide
# is therefore suppressed and panel A's key is overridden to a filled square,
# so the single collected legend reads for both panels. Colours come from
# CLUSTER_PALETTE, which is deliberately not the diagnosis palette.

FIG4_TITLE <- paste0("Figure 4. Composition-Adjusted Methylome, ",
                     "PARDS Only Cohort (n = ", n_pards, ")")
TITLE_A <- "PCA on the composition-adjusted methylome"
TITLE_B <- "Hierarchical clustering of the composition-adjusted methylome"

lab_pards       <- cluster_labels(cl_pards)
pards_group_cols <- setNames(CLUSTER_PALETTE[seq_along(lab_pards)],
                             unname(lab_pards))

# ---- Panel A: PCA ----
pca_df_fig <- as.data.frame(res_pards$pca$x) %>%
  rownames_to_column("Sample_Name") %>%
  left_join(assign_pards, by = "Sample_Name") %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID),
            by = "Sample_Name") %>%
  mutate(group_lab = unname(lab_pards[sub("^PardsGroup", "", PARDS_only_Group)]),
         label     = fig_label(Subject_ID))

panelA <- ggplot(pca_df_fig, aes(x = PC1, y = PC2, fill = group_lab, label = label)) +
  geom_point(shape = 21, size = 3, color = "grey25", stroke = 0.3) +
  geom_text_repel(size = 2.7, show.legend = FALSE, min.segment.length = 0,
                  segment.size = 0.3, segment.color = "grey55",
                  box.padding = 0.5, point.padding = 0.25,
                  max.overlaps = Inf, seed = 1) +
  scale_fill_manual(values = pards_group_cols,
                    name   = paste0("Cluster (k = ", k_compare, ")"),
                    guide  = guide_legend(
                      override.aes = list(shape = 22, size = 4.5))) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(x = paste0("PC1 (", round(res_pards$var_exp[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(res_pards$var_exp[2] * 100, 1), "%)"),
       title = TITLE_A) +
  theme_classic(base_size = 12)

# ---- Panel B: dendrogram ----
dend_b <- as.dendrogram(res_pards$hc)
leaf_b <- labels(dend_b)
grp_b  <- unname(lab_pards[as.character(cl_pards[leaf_b])])
labels(dend_b) <- fig_label(unname(setNames(sample_sheet$Subject_ID,
                                            sample_sheet$Sample_Name)[leaf_b]))

ggd_b  <- dendextend::as.ggdend(dend_b)
segs_b <- ggd_b$segments
labs_b <- ggd_b$labels %>% mutate(group_lab = grp_b)

hts_b    <- sort(res_pards$hc$height)
n_leaf_b <- length(leaf_b)
cut_h_b  <- mean(c(hts_b[n_leaf_b - k_compare], hts_b[n_leaf_b - k_compare + 1]))

max_h_b <- max(segs_b$y, na.rm = TRUE)
strip_h <- max_h_b * 0.035
strip_y <- -strip_h / 2 - max_h_b * 0.01
label_y <- -strip_h * 1.5

panelB <- ggplot() +
  geom_segment(data = segs_b, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.35, color = "black") +
  geom_hline(yintercept = cut_h_b, linetype = "dashed", color = "grey40") +
  geom_tile(data = labs_b, aes(x = x, y = strip_y, fill = group_lab),
            height = strip_h, width = 0.85) +
  geom_text(data = labs_b, aes(x = x, y = label_y, label = label),
            angle = 90, hjust = 1, vjust = 0.5, size = 2.4) +
  scale_fill_manual(values = pards_group_cols,
                    name = paste0("Cluster (k = ", k_compare, ")")) +
  guides(fill = "none") +          # panel A carries the single collected legend
  scale_x_continuous(expand = expansion(add = 0.7)) +
  scale_y_continuous(breaks = pretty(c(0, max_h_b))) +
  coord_cartesian(ylim = c(-max_h_b * 0.45, max_h_b * 1.02), clip = "off") +
  labs(x = NULL, y = "Height", title = TITLE_B) +
  theme_classic(base_size = 12) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank())

fig_panel <- (panelA / panelB) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title = FIG4_TITLE,
    theme = theme(plot.title = element_text(face = "bold", size = 13,
                                            hjust = 0.5))) &
  theme(legend.position     = "bottom",
        plot.tag            = element_text(face = "bold", size = 14),
        plot.title          = element_text(hjust = 0.5),
        plot.title.position = "plot")

ggsave(file.path(FIG_DIR, "Figure4_PARDS_only_panels.pdf"),
       fig_panel, width = 10, height = 11, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure4_PARDS_only_panels.png"),
       fig_panel, width = 10, height = 11, dpi = 600)
cat("\nFigure 4 written to:", FIG_DIR, "\n")


# -----------------------------------------------------------------------------
# IV.7 New PARDS only: adjustment re-derived within the cohort
# -----------------------------------------------------------------------------

# S5.2 Composition adjustment re-derived within the new specimens
# -----------------------------------------------------------------------------
# Same two prespecified reductions as the manuscript methods: drop the
# largest-mean cell type as reference, then drop any remaining fraction with
# variance below 1e-4.

drop_col_new <- cell_cols[which.max(colMeans(prop_mat_new))]
cat("\nReference cell type dropped (largest mean):", drop_col_new, "\n")

cand_new     <- setdiff(cell_cols, drop_col_new)
col_vars_new <- apply(prop_mat_new[, cand_new, drop = FALSE], 2, var)
nzv_new      <- names(col_vars_new)[col_vars_new < VAR_THRESHOLD]
if (length(nzv_new) > 0)
  cat("Dropped for variance <", format(VAR_THRESHOLD), ":",
      paste(nzv_new, collapse = ", "), "\n")

keep_new     <- setdiff(cand_new, nzv_new)
prop_adj_new <- prop_mat_new[, keep_new, drop = FALSE]
resid_df_new <- n_new - ncol(prop_adj_new) - 1

cat("Covariates retained:", paste(keep_new, collapse = ", "),
    "(", ncol(prop_adj_new), ")\n")
cat("Residual degrees of freedom:", resid_df_new, "\n")
if (resid_df_new < 6)
  cat("WARNING: fewer than 6 residual df. The adjustment is close to\n",
      "  saturating this subset; treat the result as descriptive only.\n")

m_adj_new <- removeBatchEffect(m_new, covariates = prop_adj_new)
ok_new    <- rowSums(is.na(m_adj_new) | !is.finite(m_adj_new)) == 0
if (sum(!ok_new) > 0)
  cat("Dropping", sum(!ok_new), "probes with non-finite residuals\n")
m_adj_new <- m_adj_new[ok_new, , drop = FALSE]
cat("Analytic matrix:", nrow(m_adj_new), "probes x", ncol(m_adj_new), "specimens\n")

# -----------------------------------------------------------------------------
# IV.8 New PARDS only: PCA, clustering, silhouette, and Table 5
# -----------------------------------------------------------------------------

# S5.3 PCA, clustering, silhouette
# -----------------------------------------------------------------------------
res_new <- cluster_diagnostics(m_adj_new, "New PARDS specimens only")

sil_three <- bind_rows(
  sil_compare,
  res_new$sil %>% mutate(analysis = paste0("New PARDS only (n = ", n_new, ")"))
) %>% dplyr::select(analysis, k, mean_silhouette)

write_csv(sil_three,
          file.path(SENS_TBL_DIR, "silhouette_by_k_three_analyses.csv"))

# Cluster sizes are reported alongside every silhouette width. A partition that
# isolates one or two specimens produces a high mean silhouette without
# indicating subgroup structure, so the width alone can mislead. Sizes are
# given for all three analytic sets, not only where they are inconvenient.

hc_by_set <- list(sample_hc_adj, res_pards$hc, res_new$hc)
names(hc_by_set) <- c(
  grep("^Full cohort", unique(sil_three$analysis), value = TRUE)[1],
  grep("^PARDS only",  unique(sil_three$analysis), value = TRUE)[1],
  grep("^New PARDS",   unique(sil_three$analysis), value = TRUE)[1])

sizes_by_k <- map_dfr(names(hc_by_set), function(nm) {
  hc <- hc_by_set[[nm]]
  if (is.null(hc)) return(NULL)
  kmax <- min(6, length(hc$order) - 1)
  map_dfr(2:kmax, function(k) {
    tb <- table(cutree(hc, k = k))
    tibble(analysis = nm, k = k,
           sizes = paste(sort(as.integer(tb), decreasing = TRUE), collapse = "/"))
  })
})

suppTable_sil3 <- sil_three %>%
  left_join(sizes_by_k, by = c("analysis", "k")) %>%
  mutate(cell = ifelse(is.na(mean_silhouette), "\u2014",
                       ifelse(is.na(sizes),
                              sprintf("%.3f", mean_silhouette),
                              sprintf("%.3f (%s)", mean_silhouette, sizes)))) %>%
  dplyr::select(analysis, k, cell) %>%
  pivot_wider(names_from = analysis, values_from = cell) %>%
  rename(`Number of clusters (k)` = k)

cat("\nSilhouette width with cluster sizes:\n")
print(as.data.frame(suppTable_sil3))

save_pub_table(
  suppTable_sil3,
  file.path(TBL_DIR, "Table5_silhouette_by_k_all_cohorts"),
  caption = paste0("Table 5. Mean silhouette width and cluster ",
                   "sizes for k = 2 to 6 in the full cohort, in PARDS ",
                   "participants only, and after excluding specimens reused ",
                   "from the previous study."),
  note = paste0("Values are the mean silhouette width with the resulting ",
                "cluster sizes in parentheses. Cell-type deconvolution, ",
                "composition adjustment, and the top 5% variance filter were ",
                "re-derived within each analytic set. The final column ",
                "excludes the ", length(DUPLICATE_SUBJECTS),
                " specimens shared with Williams et al., Respir Res 2022;23:181, ",
                "leaving ", n_new, " PARDS specimens analyzed for the first ",
                "time. Cluster sizes are reported because a partition that ",
                "isolates one or two specimens yields a high mean silhouette ",
                "width without indicating subgroup structure. Interpretation ",
                "thresholds: greater than 0.5 strong, 0.25 to 0.5 weak, less ",
                "than 0.25 essentially absent."))

# -----------------------------------------------------------------------------
# IV.9 New PARDS only: concordance, cluster sizes, per-sample silhouette
# -----------------------------------------------------------------------------

# S5.4 Concordance with the PARDS-only partition from S1
# -----------------------------------------------------------------------------
# Do the new specimens group the same way when the reused specimens are gone?

k_new <- min(k_compare, n_new - 1L)
cl_new <- cutree(res_new$hc, k = k_new)
stopifnot(identical(names(cl_new), labels(res_new$dist)))

conc_new <- tibble(
  Sample_Name   = names(cl_new),
  New_Only_Group = paste0("NewGroup", unname(cl_new))
) %>%
  left_join(assign_pards, by = "Sample_Name") %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID),
            by = "Sample_Name") %>%
  dplyr::select(Subject_ID, PARDS_only_Group, New_Only_Group) %>%
  arrange(Subject_ID)

cat("\nAssignment at k =", k_new, ": PARDS-only (S1) vs new-specimens-only\n")
print(table(conc_new$PARDS_only_Group, conc_new$New_Only_Group,
            dnn = c("PARDS only (n = 21)", paste0("New only (n = ", n_new, ")"))))

ari_new <- adjusted_rand_index(conc_new$PARDS_only_Group, conc_new$New_Only_Group)
cat("Adjusted Rand index vs the S1 PARDS-only partition:", round(ari_new, 3), "\n")

write_csv(conc_new,
          file.path(SENS_TBL_DIR, "cluster_concordance_new_specimens_only.csv"))
# S5.4b Cluster sizes and per-sample silhouette
# -----------------------------------------------------------------------------
# A silhouette value is uninterpretable without the cluster sizes beside it: a
# lone specimen split off from a compact remainder inflates the mean without
# any subgroup being present. Both are reported here and carried into the
# figure legend.

cat("\nCluster sizes at k =", k_new, ":\n"); print(table(cl_new))

size_by_k_new <- map_dfr(2:min(6, n_new - 1), function(k) {
  tb <- table(cutree(res_new$hc, k = k))
  tibble(k = k,
         cluster_sizes = paste(sort(as.integer(tb), decreasing = TRUE),
                               collapse = " / "),
         smallest = min(as.integer(tb)))
})
cat("\nCluster sizes by k (new specimens only):\n")
print(as.data.frame(size_by_k_new))
write_csv(size_by_k_new,
          file.path(SENS_TBL_DIR, "cluster_sizes_new_specimens_only.csv"))

sil_new_obj <- silhouette(cl_new, res_new$dist)
sil_new_tbl <- tibble(
  Sample_Name = names(cl_new),
  cluster     = paste0("Group ", sil_new_obj[, 1]),
  sil_width   = round(sil_new_obj[, 3], 3)
) %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID),
            by = "Sample_Name") %>%
  dplyr::select(Subject_ID, cluster, sil_width) %>%
  arrange(cluster, desc(sil_width))

cat("\nPer-sample silhouette (new specimens only, k =", k_new, "):\n")
print(as.data.frame(sil_new_tbl))
write_csv(sil_new_tbl,
          file.path(SENS_TBL_DIR, "per_sample_silhouette_new_specimens_only.csv"))

singleton_k <- names(table(cl_new))[table(cl_new) == 1]
if (length(singleton_k) > 0) {
  singleton_id <- sample_sheet$Subject_ID[
    match(names(cl_new)[cl_new %in% as.integer(singleton_k)],
          sample_sheet$Sample_Name)]
  cat("\nSingleton cluster(s) at k =", k_new, ":",
      paste(singleton_id, collapse = ", "), "\n")
  cat("Epithelial fraction of singleton(s):",
      paste(sprintf("%.3f", prop_mat_new[
        names(cl_new)[cl_new %in% as.integer(singleton_k)], "Epi"]),
        collapse = ", "), "\n")
}

# -----------------------------------------------------------------------------
# IV.10 New PARDS only: figure (Figure 5)
# -----------------------------------------------------------------------------

# S5.5b Paneled figure matching the PARDS-only figure
# -----------------------------------------------------------------------------
# Same layout and palette as the PARDS-only panels: PCA above, dendrogram
# below, one shared legend. Cluster sizes are written into the legend labels so
# the singleton is visible without reading the text.

FIG5_TITLE <- paste0("Figure 5. Composition-Adjusted Methylome, New PARDS ",
                     "Only Cohort, Excluding Reused Specimens (n = ", n_new, ")")
TITLE_A_NEW <- "PCA on the composition-adjusted methylome"
TITLE_B_NEW <- "Hierarchical clustering of the composition-adjusted methylome"

lab_new_base   <- cluster_labels(cl_new)
lab_new        <- setNames(lab_new_base,
                           paste0("NewGroup", names(lab_new_base)))
new_group_cols <- setNames(CLUSTER_PALETTE[seq_along(lab_new)],
                           unname(lab_new))

# ---- Panel A: PCA ----
pcaA_new <- as.data.frame(res_new$pca$x) %>%
  rownames_to_column("Sample_Name") %>%
  mutate(grp_key   = paste0("NewGroup", cl_new[Sample_Name]),
         group_lab = unname(lab_new[grp_key])) %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID),
            by = "Sample_Name") %>%
  mutate(label = fig_label(Subject_ID))

panelA_new <- ggplot(pcaA_new,
                     aes(x = PC1, y = PC2, fill = group_lab, label = label)) +
  geom_point(shape = 21, size = 3, color = "grey25", stroke = 0.3) +
  geom_text_repel(size = 2.7, show.legend = FALSE, min.segment.length = 0,
                  segment.size = 0.3, segment.color = "grey55",
                  box.padding = 0.5, point.padding = 0.25,
                  max.overlaps = Inf, seed = 1) +
  scale_fill_manual(values = new_group_cols,
                    name   = paste0("Cluster (k = ", k_new, ")"),
                    guide  = guide_legend(
                      override.aes = list(shape = 22, size = 4.5))) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(x = paste0("PC1 (", round(res_new$var_exp[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(res_new$var_exp[2] * 100, 1), "%)"),
       title = TITLE_A_NEW) +
  theme_classic(base_size = 12)

# ---- Panel B: dendrogram ----
dend_n <- as.dendrogram(res_new$hc)
leaf_n <- labels(dend_n)
grp_n  <- unname(lab_new[paste0("NewGroup", cl_new[leaf_n])])
labels(dend_n) <- fig_label(unname(setNames(sample_sheet$Subject_ID,
                                            sample_sheet$Sample_Name)[leaf_n]))

ggd_n  <- dendextend::as.ggdend(dend_n)
segs_n <- ggd_n$segments
labs_n <- ggd_n$labels %>% mutate(group_lab = grp_n)

hts_n   <- sort(res_new$hc$height)
nleaf_n <- length(leaf_n)
cut_h_n <- mean(c(hts_n[nleaf_n - k_new], hts_n[nleaf_n - k_new + 1]))

max_h_n  <- max(segs_n$y, na.rm = TRUE)
strip_hn <- max_h_n * 0.035
strip_yn <- -strip_hn / 2 - max_h_n * 0.01
label_yn <- -strip_hn * 1.5

panelB_new <- ggplot() +
  geom_segment(data = segs_n, aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.35, color = "black") +
  geom_hline(yintercept = cut_h_n, linetype = "dashed", color = "grey40") +
  geom_tile(data = labs_n, aes(x = x, y = strip_yn, fill = group_lab),
            height = strip_hn, width = 0.85) +
  geom_text(data = labs_n, aes(x = x, y = label_yn, label = label),
            angle = 90, hjust = 1, vjust = 0.5, size = 2.4) +
  scale_fill_manual(values = new_group_cols,
                    name = paste0("Cluster (k = ", k_new, ")")) +
  guides(fill = "none") +          # panel A carries the single collected legend
  scale_x_continuous(expand = expansion(add = 0.7)) +
  scale_y_continuous(breaks = pretty(c(0, max_h_n))) +
  coord_cartesian(ylim = c(-max_h_n * 0.45, max_h_n * 1.02), clip = "off") +
  labs(x = NULL, y = "Height", title = TITLE_B_NEW) +
  theme_classic(base_size = 12) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x  = element_blank())

fig_panel_new <- (panelA_new / panelB_new) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title = FIG5_TITLE,
    theme = theme(plot.title = element_text(face = "bold", size = 13,
                                            hjust = 0.5))) &
  theme(legend.position     = "bottom",
        plot.tag            = element_text(face = "bold", size = 14),
        plot.title          = element_text(hjust = 0.5),
        plot.title.position = "plot")

ggsave(file.path(FIG_DIR, "Figure5_new_specimens_only_panels.pdf"),
       fig_panel_new, width = 10, height = 11, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure5_new_specimens_only_panels.png"),
       fig_panel_new, width = 10, height = 11, dpi = 600)
cat("\nFigure 5 written to:", FIG_DIR, "\n")

# -----------------------------------------------------------------------------
# IV.11 Supplemental Table 3: PCA variance before and after adjustment
# -----------------------------------------------------------------------------

# composition adjustment, in each cohort. Both are computed here rather than
# assembled from S9 and S6.4, so the two halves of the comparison come from one
# code path and cannot drift apart.
#
# The adjustment is re-derived within each cohort, following the rule stated in
# the Methods: the cell type with the largest mean fraction is dropped as the
# reference, and any remaining fraction with variance below 1e-4 is dropped.

cat("\n=== PART IV.11. PC1 before and after composition adjustment ===\n")

pc_top2 <- function(m) {
  pv    <- rowVars(m, na.rm = TRUE)
  m_top <- m[pv >= quantile(pv, 0.95, na.rm = TRUE), , drop = FALSE]
  pca   <- prcomp(t(m_top), center = TRUE, scale. = FALSE)
  ve    <- (pca$sdev^2) / sum(pca$sdev^2)
  c(pc1 = ve[1] * 100, pc2 = ve[2] * 100, probes = nrow(m_top))
}

pc_pre_post <- map_dfr(ANALYTIC_SETS, function(s) {
  pre <- pc_top2(s$m)

  ref  <- cell_cols[which.max(colMeans(s$prop))]
  cand <- setdiff(cell_cols, ref)
  v    <- apply(s$prop[, cand, drop = FALSE], 2, var)
  keep <- cand[v >= VAR_THRESHOLD]
  cov  <- s$prop[, keep, drop = FALSE]

  m_adj <- removeBatchEffect(s$m, covariates = cov)
  ok    <- rowSums(is.na(m_adj) | !is.finite(m_adj)) == 0
  post  <- pc_top2(m_adj[ok, , drop = FALSE])

  tibble(cohort        = s$label,
         n             = ncol(s$m),
         n_covariates  = ncol(cov),
         residual_df   = ncol(s$m) - ncol(cov) - 1,
         pc1_pre       = unname(pre["pc1"]),
         pc1_post      = unname(post["pc1"]),
         pc1_reduction = unname(pre["pc1"] - post["pc1"]),
         pc2_pre       = unname(pre["pc2"]),
         pc2_post      = unname(post["pc2"]))
})

cat("\nPC1 before and after adjustment:\n")
print(as.data.frame(pc_pre_post %>% mutate(across(where(is.numeric), ~round(.x, 1)))))
write_csv(pc_pre_post,
          file.path(TBL_DIR, "SuppTable3_pca_variance_pre_post_raw.csv"))

save_pub_table(
  pc_pre_post %>%
    transmute(`Cohort`                        = cohort,
              `Specimens (n)`                 = n,
              `Covariates (n)`                = n_covariates,
              `Residual df`                   = residual_df,
              `PC1 before adjustment (%)`     = sprintf("%.1f", pc1_pre),
              `PC1 after adjustment (%)`      = sprintf("%.1f", pc1_post),
              `PC2 before adjustment (%)`     = sprintf("%.1f", pc2_pre),
              `PC2 after adjustment (%)`      = sprintf("%.1f", pc2_post)),
  file.path(TBL_DIR, "SuppTable3_pca_variance_pre_post"),
  caption = paste0("Supplemental Table 3. Variance explained by the first two ",
                   "principal components before and after adjustment for ",
                   "cellular composition."),
  note = paste0("Principal component analysis was performed on the top 5% ",
                "most variable CpGs, reselected within each cohort and ",
                "separately before and after adjustment. The composition ",
                "adjustment was re-derived independently within each cohort ",
                "using the rule stated in the Methods: the cell type with the ",
                "largest mean fraction was dropped as the reference and any ",
                "remaining fraction with variance below 1e-4 was dropped. ",
                "Residual degrees of freedom are the specimens remaining after ",
                "the covariates and intercept are fitted, and fall to ",
                min(pc_pre_post$residual_df),
                " in the smallest cohort."))

saveRDS(pc_pre_post, file.path(RDS_DIR, "checkpoint_pc1_pre_post_S9b.rds"))

cat("\n--- For the Results text ---\n")
for (i in seq_len(nrow(pc_pre_post))) {
  r <- pc_pre_post[i, ]
  cat(sprintf("%s: PC1 fell from %.1f%% to %.1f%% (reduction of %.1f points); PC2 %.1f%% to %.1f%%.\n",
              r$cohort, r$pc1_pre, r$pc1_post, r$pc1_reduction, r$pc2_pre, r$pc2_post))
}

cat("\nPART IV complete:", format(Sys.time()), "\n")


# =============================================================================


# =============================================================================
# PART V. CELL-TYPE-SPECIFIC DIFFERENTIAL METHYLATION (CellDMC)
# =============================================================================

# Two goals:
#   S3.1 Determine definitively whether any cell-type compartment was dropped
#        by the near-invariance rule in the immunocompromised contrast, so the
#        methods can state it rather than describe a rule of unknown effect.
#   S3.2 Re-run the permutation null at B = 1000 (master section 15c used 100,
#        which floors the empirical p at 1/101 = 0.0099).
#
# Controls are retained in this contrast by design. The subgroup-discovery
# objections from Reviewers 1 and 3 concern clustering within PARDS; the
# immunocompromised contrast is a prespecified clinical comparison across the
# full cohort, and the PARDS-only sensitivity analysis in S1 showed removing
# controls did not change the compositional structure.
#
# Self-contained: rebuilds frac_cdmc and the IC phenotype from the checkpoint
# objects rather than depending on section 15b still being in the session.

cat("\n=== PART V. CellDMC in the immunocompromised contrast ===\n")

suppressPackageStartupMessages(library(EpiDISH))

# -----------------------------------------------------------------------------
# S3.1 Rebuild the CellDMC fraction matrix and check the invariance rule
# -----------------------------------------------------------------------------
# Same 3-compartment collapse as master section 15b: Epi, Neutro, and Rest as
# the complement of the two. Because Rest is defined as 1 - Epi - Neutro, the
# three sum to one by construction and the rowSums division below only guards
# floating-point drift; it is not a substantive renormalization step.

frac_full_s3 <- as.matrix(
  cell_props[match(colnames(beta_collapsed), cell_props$Sample_Name),
             cell_cols, drop = FALSE]
)
rownames(frac_full_s3) <- colnames(beta_collapsed)

Epi_s3    <- frac_full_s3[, "Epi"]
Neutro_s3 <- frac_full_s3[, "Neutro"]
Rest_s3   <- pmax(0, 1 - Epi_s3 - Neutro_s3)
frac_cdmc_s3 <- cbind(Epi = Epi_s3, Neutro = Neutro_s3, Rest = Rest_s3)

drift <- max(abs(rowSums(frac_cdmc_s3) - 1))
cat("Max deviation of compartment row sums from 1 before guard:",
    signif(drift, 3), "\n")
frac_cdmc_s3 <- frac_cdmc_s3 / rowSums(frac_cdmc_s3)

# IC phenotype, aligned to the columns of beta_collapsed.
# Immunocompromised_derived is created in master section 13, which runs AFTER
# the section-12 checkpoint this script resumes from, so the column is absent
# when running standalone. Derive it here using the identical rule.
if (!"Immunocompromised_derived" %in% names(sample_sheet)) {
  if (!"Principal_Comorbidity" %in% names(sample_sheet))
    stop("Neither Immunocompromised_derived nor Principal_Comorbidity is present ",
         "in sample_sheet; cannot define the immunocompromised contrast.")
  cat("Deriving Immunocompromised_derived from Principal_Comorbidity ",
      "(same rule as master section 13).\n", sep = "")
  sample_sheet <- sample_sheet %>%
    mutate(Immunocompromised_derived = grepl("immun|onc|transplant|hsct|chemo",
                                             Principal_Comorbidity,
                                             ignore.case = TRUE))
}

ss_s3 <- sample_sheet[match(rownames(frac_cdmc_s3), sample_sheet$Sample_Name), ]
stopifnot("Immunocompromised_derived" %in% names(ss_s3))
ic_pheno_s3 <- as.integer(ss_s3$Immunocompromised_derived)
names(ic_pheno_s3) <- rownames(frac_cdmc_s3)
stopifnot(!any(is.na(ic_pheno_s3)))

# The manuscript defines immunocompromised as current chemotherapy, status post
# bone marrow transplant, or an active immunosuppressive regimen. The regex
# above operationalizes that from the recorded principal comorbidity. Print the
# mapping so the coded definition can be checked against the stated one.
cat("\nPrincipal comorbidity by derived immunocompromised flag:\n")
print(ss_s3 %>%
        dplyr::select(Subject_ID, Principal_Comorbidity,
                      Immunocompromised_derived) %>%
        arrange(desc(Immunocompromised_derived), Subject_ID) %>%
        as.data.frame())

cat("IC contrast group sizes: immunocompromised n =", sum(ic_pheno_s3 == 1L),
    ", non-immunocompromised n =", sum(ic_pheno_s3 == 0L), "\n")
cat("Controls included in this contrast:",
    sum(ss_s3$Sample_Group == CONTROL_LABEL), "\n")
cat("Immunocompromised subjects that are controls:",
    sum(ic_pheno_s3 == 1L & ss_s3$Sample_Group == CONTROL_LABEL), "\n")

# The invariance rule, evaluated at BOTH the threshold coded in master section
# 15b (1e-5) and the threshold stated in the manuscript methods (1e-4). These
# differ by an order of magnitude and must be reconciled before submission.
THRESH_CODED  <- 1e-5
THRESH_METHODS <- 1e-4

v_compartment <- apply(frac_cdmc_s3, 2, var)

invariance_check <- tibble(
  compartment       = names(v_compartment),
  mean_prop         = colMeans(frac_cdmc_s3),
  variance          = as.numeric(v_compartment),
  n_nonzero         = colSums(frac_cdmc_s3 > 1e-8),
  dropped_at_1e5    = v_compartment < THRESH_CODED,
  dropped_at_1e4    = v_compartment < THRESH_METHODS
) %>% arrange(desc(mean_prop))

cat("\nNear-invariance check, immunocompromised contrast:\n")
print(invariance_check %>%
        mutate(mean_pct = signif(mean_prop * 100, 3),
               variance = signif(variance, 3)) %>%
        dplyr::select(compartment, mean_pct, variance, n_nonzero,
                      dropped_at_1e5, dropped_at_1e4) %>%
        as.data.frame())

n_drop_coded   <- sum(invariance_check$dropped_at_1e5)
n_drop_methods <- sum(invariance_check$dropped_at_1e4)

cat("\nCompartments dropped at the coded threshold (", THRESH_CODED, "): ",
    n_drop_coded, "\n", sep = "")
cat("Compartments dropped at the methods threshold (", THRESH_METHODS, "): ",
    n_drop_methods, "\n", sep = "")
if (n_drop_coded == 0 && n_drop_methods == 0) {
  cat("=> No compartment was dropped under either threshold. The methods can\n",
      "   state that the rule was applied and did not remove any compartment,\n",
      "   and the choice of threshold is immaterial to the reported result.\n")
} else {
  cat("=> At least one compartment is affected. The threshold discrepancy\n",
      "   between the code and the methods text is material and must be\n",
      "   resolved before the results are reported.\n")
}

write_csv(invariance_check,
          file.path(SENS_TBL_DIR, "celldmc_IC_invariance_check.csv"))


# -----------------------------------------------------------------------------
# V.2 Observed CellDMC fit
# -----------------------------------------------------------------------------
# One genome-wide fit, roughly 8 minutes.

cat("\nComputing the observed CellDMC result (unadjusted, all compartments)...\n")
res_obs_s3 <- CellDMC(beta.m = beta_collapsed, pheno.v = ic_pheno_s3,
                      frac.m = frac_cdmc_s3, adjPMethod = "fdr",
                      adjPThresh = 0.05, mc.cores = 1)
obs_counts_s3 <- colSums(abs(res_obs_s3$dmct[, -1, drop = FALSE]) > 0,
                         na.rm = TRUE)
cat("Observed DMCs (FDR < 0.05):\n"); print(obs_counts_s3)
saveRDS(obs_counts_s3, file.path(RDS_DIR, "celldmc_IC_observed_counts.rds"))

# -----------------------------------------------------------------------------
# V.3 Permutation null, B = 1000
# -----------------------------------------------------------------------------
# RUNTIME. This reruns CellDMC genome-wide B times. At roughly 8 minutes per
# fit, B = 1000 is not feasible on a laptop; the reported null was computed on
# the UAMS Grace cluster with a 20-task SLURM array. The scripts that did it
# are in hpc/: celldmc_perm_array.sbatch submits the array,
# celldmc_perm_array.R computes one slice, and celldmc_perm_combine.R merges
# the slices and writes the summary. The merged result is copied back to
# RDS_DIR, where the branch below picks it up.
#
# Behaviour:
#   - if the cluster result is present in RDS_DIR, it is loaded and used;
#   - otherwise, if RUN_PERMUTATION is TRUE, the null is computed here, with a
#     checkpoint written every CHUNK permutations so an interrupted run resumes
#     rather than restarts;
#   - otherwise the section is skipped and Table 6 reports the null as missing.

RUN_PERMUTATION <- TRUE
PERM_CORES      <- max(1L, parallel::detectCores() - 2L)
CHUNK           <- 50

# -----------------------------------------------------------------------------
# Everything the cluster array needs, written where hpc/celldmc_perm_array.R
# expects it. Writing this here rather than assembling it by hand is what makes
# the cluster path reproducible: the observed counts travel with the data, so
# the permutation null can never be anchored to counts from a different run.
# -----------------------------------------------------------------------------
perm_input_file <- file.path(RDS_DIR, "celldmc_perm_input.rds")
saveRDS(list(beta_collapsed = beta_collapsed,
             frac_cdmc      = frac_cdmc_s3,
             ic_pheno       = ic_pheno_s3,
             obs_counts     = obs_counts_s3,
             fdr_thresh     = 0.05,
             b_perm         = B_PERM,
             seed           = 20260101,
             rng_kind       = "L'Ecuyer-CMRG",
             created        = format(Sys.time())),
        perm_input_file)
cat("Cluster input written to:", perm_input_file, "\n")
cat("  To run the null on a cluster, copy it next to the scripts in hpc/ and\n")
cat("  submit hpc/celldmc_perm_array.sbatch; see hpc/README section in README.md.\n")

# The cluster array uses L'Ecuyer-CMRG so that forked workers get independent
# substreams. Match it here, so permutation b is the same relabelling whether it
# was computed locally or on the cluster, and slices from the two can be pooled.
RNGkind("L'Ecuyer-CMRG")

perm_cluster_file <- file.path(RDS_DIR,
                               "celldmc_IC_permutation_null_B1000_full.rds")
perm_file <- file.path(RDS_DIR, "celldmc_IC_permutation_null_B1000.rds")

perm_summary_s3 <- NULL

if (file.exists(perm_cluster_file)) {

  perm_cluster    <- readRDS(perm_cluster_file)
  perm_summary_s3 <- as_tibble(perm_cluster$perm_summary)
  perm_mat_s3     <- perm_cluster$perm_mat
  perm_ok_s3      <- perm_mat_s3[complete.cases(perm_mat_s3), , drop = FALSE]
  min_p           <- perm_cluster$min_p
  stopifnot(all(perm_summary_s3$observed ==
                  obs_counts_s3[perm_summary_s3$cell_type]))
  cat("\nPermutation null loaded from the cluster run: ", nrow(perm_ok_s3),
      " of ", nrow(perm_mat_s3), " permutations completed\n", sep = "")

} else if (RUN_PERMUTATION) {

  cat("\nCluster result not found; computing the permutation null locally.\n")
  cat("Expect roughly", round(B_PERM * 8 / max(1L, PERM_CORES) / 60),
      "hours on", PERM_CORES, "cores.\n")

  perm_state <- if (file.exists(perm_file)) readRDS(perm_file) else NULL
  if (!is.null(perm_state) &&
      !identical(dim(perm_state$perm_mat),
                 c(B_PERM, ncol(frac_cdmc_s3)))) {
    cat("Existing checkpoint has different dimensions; starting fresh.\n")
    perm_state <- NULL
  }
  if (is.null(perm_state)) {
    perm_state <- list(
      perm_mat  = matrix(NA_integer_, nrow = B_PERM, ncol = ncol(frac_cdmc_s3),
                         dimnames = list(NULL, colnames(frac_cdmc_s3))),
      last_done = 0L, seed = 20260101)
  } else {
    cat("Resuming from checkpoint:", perm_state$last_done, "of", B_PERM, "\n")
  }

  set.seed(perm_state$seed)
  # Draw and discard the permutations already completed so a resumed run
  # continues the same random sequence rather than repeating it.
  if (perm_state$last_done > 0)
    for (b in seq_len(perm_state$last_done)) invisible(sample(ic_pheno_s3))

  t_start    <- Sys.time()
  start_from <- perm_state$last_done + 1L
  if (perm_state$last_done < B_PERM) {
    for (b in seq.int(start_from, B_PERM)) {
      perm_pheno <- sample(ic_pheno_s3)        # preserves the group sizes
      rp <- tryCatch(
        CellDMC(beta.m = beta_collapsed, pheno.v = perm_pheno,
                frac.m = frac_cdmc_s3, adjPMethod = "fdr",
                adjPThresh = 0.05, mc.cores = PERM_CORES),
        error = function(e) NULL)
      if (!is.null(rp))
        perm_state$perm_mat[b, ] <- colSums(abs(rp$dmct[, -1, drop = FALSE]) > 0,
                                            na.rm = TRUE)
      perm_state$last_done <- b
      if (b %% CHUNK == 0 || b == B_PERM) {
        saveRDS(perm_state, perm_file)
        elapsed_min  <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
        min_per_perm <- elapsed_min / (b - start_from + 1L)
        cat(sprintf("  permutation %d / %d | %.1f min elapsed | %.2f min each | ~%.0f min remaining\n",
                    b, B_PERM, elapsed_min, min_per_perm,
                    min_per_perm * (B_PERM - b)))
      }
    }
  }

  perm_mat_s3 <- perm_state$perm_mat
  perm_ok_s3  <- perm_mat_s3[complete.cases(perm_mat_s3), , drop = FALSE]
  min_p       <- 1 / (nrow(perm_ok_s3) + 1)
  cat("\nPermutations completed:", nrow(perm_ok_s3), "of", B_PERM, "\n")

  perm_summary_s3 <- map_dfr(colnames(perm_mat_s3), function(ct) {
    pc <- perm_ok_s3[, ct]
    tibble(cell_type   = ct,
           observed    = obs_counts_s3[ct],
           perm_median = median(pc),
           perm_mean   = mean(pc),
           perm_p95    = unname(quantile(pc, 0.95)),
           perm_max    = max(pc),
           # add-one smoothing so the empirical p is never exactly zero
           emp_p_value = (sum(pc >= obs_counts_s3[ct]) + 1) / (length(pc) + 1))
  })

} else {
  cat("\nPermutation null skipped (RUN_PERMUTATION = FALSE and no cluster",
      "result found).\n")
}

if (!is.null(perm_summary_s3)) {
  cat("\nPermutation null summary (B =", nrow(perm_ok_s3), "):\n")
  print(perm_summary_s3 %>% mutate(across(where(is.numeric), ~round(.x, 4))))
  cat("Smallest achievable empirical p at this B:", signif(min_p, 3), "\n")
  write_csv(perm_summary_s3,
            file.path(SENS_TBL_DIR, "celldmc_IC_permutation_null_B1000.csv"))
}


# -----------------------------------------------------------------------------
# V.4 Covariate-adjusted rerun (age + slide)
# -----------------------------------------------------------------------------

# including one immunocompromised case, so including it silently reduced the
# contrast below the stated 5 vs 19. With age and slide only, all 24 specimens
# should be retained; confirm the printed counts are 5 and 19.
#
# One genome-wide CellDMC fit, roughly 8 minutes.

cat("\n--- S3.4 Covariate-adjusted rerun (age + slide) ---\n")

# The slide/array identifier is not named consistently across sample sheets,
# and it is not present in the section-12 checkpoint under "Slide". Resolve it
# from candidates; if none is present, fall back to an age-only model and say
# so, rather than silently dropping specimens or failing.
SLIDE_CANDIDATES <- c("Slide", "Sentrix_ID", "SentrixID", "Sentrix_Barcode",
                      "Slide_ID", "Chip", "BeadChip", "Beadchip", "Array",
                      "Array_ID", "Barcode", "Basename")
slide_col <- SLIDE_CANDIDATES[SLIDE_CANDIDATES %in% names(ss_s3)][1]

if (is.na(slide_col)) {
  cat("No slide/array column found among:",
      paste(SLIDE_CANDIDATES, collapse = ", "), "\n")
  cat("Columns present in sample_sheet:\n"); print(names(ss_s3))
  cat("Falling back to an age-only adjusted model. Add the correct name to\n",
      "SLIDE_CANDIDATES to include slide, and update the methods text if the\n",
      "adjusted model is reported without it.\n")
  slide_vals <- NULL
} else {
  cat("Using slide/array column:", slide_col, "\n")
  slide_vals <- factor(ss_s3[[slide_col]])
  cat("  distinct slides:", nlevels(slide_vals), "\n")
  if (nlevels(slide_vals) >= ncol(beta_collapsed) - 3L)
    cat("  WARNING: nearly one slide per specimen; this term will consume\n",
        "  most of the residual degrees of freedom.\n")
}

cov_df_s3 <- if (is.null(slide_vals)) {
  data.frame(Age = suppressWarnings(as.numeric(ss_s3$Age)))
} else {
  data.frame(Age = suppressWarnings(as.numeric(ss_s3$Age)), Slide = slide_vals)
}
cc_s3 <- complete.cases(cov_df_s3)

cat("Complete-covariate specimens:", sum(cc_s3), "of", nrow(cov_df_s3), "\n")
cat("  immunocompromised retained:", sum(ic_pheno_s3[cc_s3] == 1L),
    "| non-immunocompromised retained:", sum(ic_pheno_s3[cc_s3] == 0L), "\n")
if (sum(cc_s3) < nrow(cov_df_s3))
  cat("  NOTE: specimens dropped for missing covariates:",
      paste(ss_s3$Subject_ID[!cc_s3], collapse = ", "), "\n")

cov_form_s3 <- if (is.null(slide_vals)) ~ Age else ~ Age + Slide
cov_mod_s3  <- model.matrix(cov_form_s3,
                            data = cov_df_s3[cc_s3, , drop = FALSE])[, -1, drop = FALSE]
cat("Adjusted model:", deparse(cov_form_s3),
    "|", ncol(cov_mod_s3), "covariate columns |",
    sum(cc_s3) - ncol(cov_mod_s3) - 2L, "residual df\n")

res_adj_s3 <- CellDMC(beta.m  = beta_collapsed[, cc_s3, drop = FALSE],
                      pheno.v = ic_pheno_s3[cc_s3],
                      frac.m  = frac_cdmc_s3[cc_s3, , drop = FALSE],
                      cov.mod = cov_mod_s3,
                      adjPMethod = "fdr", adjPThresh = 0.05, mc.cores = 1)
adj_counts_s3 <- colSums(abs(res_adj_s3$dmct[, -1, drop = FALSE]) > 0, na.rm = TRUE)

cat("\nCovariate-adjusted DMCs (FDR < 0.05):\n"); print(adj_counts_s3)
cat("Unadjusted, for comparison:\n"); print(obs_counts_s3)


# -----------------------------------------------------------------------------
# V.5 Leave-one-out over the immunocompromised subjects
# -----------------------------------------------------------------------------

# Each immunocompromised subject is removed in turn, giving five re-analyses of
# 4 vs 19 subjects. Reviewer 2 comment 5 reads this as a 4 vs 1 comparison, so
# the group sizes are printed and carried into the table explicitly.
#
# Five genome-wide fits, run in parallel. Set LOO_CORES to 1 if memory is tight.

cat("\n--- S3.5 Leave-one-out over immunocompromised subjects ---\n")

LOO_CORES <- min(5L, max(1L, parallel::detectCores() - 2L))

ic_idx  <- which(ic_pheno_s3 == 1L)
ic_subj <- ss_s3$Subject_ID[ic_idx]
cat("Immunocompromised subjects:", paste(ic_subj, collapse = ", "), "\n")
cat("Each fold compares", length(ic_idx) - 1L, "immunocompromised vs",
    sum(ic_pheno_s3 == 0L), "non-immunocompromised subjects\n")
cat("Running", length(ic_idx), "fits on", LOO_CORES, "cores...\n")

loo_list <- parallel::mclapply(seq_along(ic_idx), function(j) {
  keep <- setdiff(seq_along(ic_pheno_s3), ic_idx[j])
  r <- tryCatch(
    CellDMC(beta.m = beta_collapsed[, keep, drop = FALSE],
            pheno.v = ic_pheno_s3[keep],
            frac.m = frac_cdmc_s3[keep, , drop = FALSE],
            adjPMethod = "fdr", adjPThresh = 0.05, mc.cores = 1),
    error = function(e) NULL)
  if (is.null(r)) return(NULL)
  colSums(abs(r$dmct[, -1, drop = FALSE]) > 0, na.rm = TRUE)
}, mc.cores = LOO_CORES)

loo_tbl <- map_dfr(seq_along(ic_idx), function(j) {
  cnt <- loo_list[[j]]
  if (is.null(cnt))
    cnt <- setNames(rep(NA_integer_, ncol(frac_cdmc_s3)), colnames(frac_cdmc_s3))
  tibble(subject_removed = ic_subj[j],
         n_immunocompromised = length(ic_idx) - 1L,
         n_non_immunocompromised = sum(ic_pheno_s3 == 0L),
         !!!as.list(cnt))
})

cat("\nLeave-one-out DMC counts:\n"); print(as.data.frame(loo_tbl))
cat("\nFull-cohort counts for comparison:\n"); print(obs_counts_s3)

write_csv(loo_tbl, file.path(SENS_TBL_DIR, "celldmc_IC_leave_one_out.csv"))

save_support_table(
  loo_tbl %>%
    rename(`Subject removed`           = subject_removed,
           `Immunocompromised (n)`     = n_immunocompromised,
           `Non-immunocompromised (n)` = n_non_immunocompromised),
  file.path(SENS_TBL_DIR, "Supporting_CellDMC_leave_one_out"),
  caption = paste0("Supporting table. Leave-one-out analysis of the ",
                   "immunocompromised contrast in CellDMC."),
  note = paste0("Each immunocompromised subject was removed in turn, yielding ",
                length(ic_idx), " re-analyses of ", length(ic_idx) - 1L,
                " immunocompromised versus ", sum(ic_pheno_s3 == 0L),
                " non-immunocompromised subjects. Values are counts of ",
                "differentially methylated cytosines at a Benjamini-Hochberg ",
                "false discovery rate below 0.05. Full-cohort counts were ",
                paste(sprintf("%s = %d", names(obs_counts_s3), obs_counts_s3),
                      collapse = ", "), "."))


# -----------------------------------------------------------------------------
# V.6 Figure 6: leave-one-out DMC counts
# -----------------------------------------------------------------------------

# Built from loo_tbl so the figure and the summary table come from one source.
# Counts are plotted on a log scale with +1 offset because the range spans two
# orders of magnitude. Dashed lines mark the full-cohort counts.
#
# Note the direction: removing a subject INFLATES the counts, so that subject's
# presence suppresses them. The driving subject differs by compartment, so it
# is identified programmatically rather than named in a fixed caption.

# Shared with Table 6 so the figure and the table cannot name the same
# compartments differently (they previously read "Rest (pooled)" vs "Others").
comp_labels <- COMPARTMENT_LABELS
comps <- intersect(names(comp_labels), names(loo_tbl))

loo_long <- loo_tbl %>%
  dplyr::select(subject_removed, all_of(comps)) %>%
  pivot_longer(all_of(comps), names_to = "compartment", values_to = "dmc") %>%
  mutate(compartment = factor(comp_labels[compartment],
                              levels = unname(comp_labels[comps])),
         subject_label = fig_label(subject_removed))

obs_line <- tibble(
  compartment = factor(unname(comp_labels[comps]), levels = unname(comp_labels[comps])),
  observed    = as.numeric(obs_counts_s3[comps]))

# Which dropped subject maximally inflates each compartment
drivers <- loo_long %>%
  group_by(compartment) %>%
  slice_max(dmc, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(obs_line, by = "compartment") %>%
  mutate(fold = dmc / pmax(observed, 1))

cat("\nLargest leave-one-out inflation by compartment:\n")
print(as.data.frame(drivers %>%
                      transmute(compartment, subject_removed, dmc, observed,
                                fold = round(fold, 1))))

driver_txt <- paste(
  sprintf("%s: %s (%.0f-fold)", drivers$compartment,
          fig_label(drivers$subject_removed), drivers$fold),
  collapse = "; ")

p_loo <- ggplot(loo_long, aes(x = subject_label, y = dmc + 1)) +
  geom_hline(data = obs_line, aes(yintercept = observed + 1),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_col(width = 0.65, fill = FIG1_BLUE, color = "grey25", linewidth = 0.3) +
  facet_wrap(~ compartment, nrow = 1) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000),
                labels = c("1", "10", "100", "1,000", "10,000"),
                expand = expansion(mult = c(0, 0.08))) +
  annotation_logticks(sides = "l", size = 0.3) +
  labs(x = "Subject removed", y = "DMC count + 1 (log scale)",
       title = paste0("Figure 6. Leave-One-Out Analysis of the ",
                      "Immunocompromised Contrast"),
       subtitle = paste0("Dashed lines mark full-cohort counts. Removing one ",
                         "subject inflates counts by up to two orders of ",
                         "magnitude.")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
        plot.title.position = "plot")

ggsave(file.path(FIG_DIR, "Figure6_CellDMC_leave_one_out.pdf"),
       p_loo, width = 9, height = 5, device = cairo_pdf)
ggsave(file.path(FIG_DIR, "Figure6_CellDMC_leave_one_out.png"),
       p_loo, width = 9, height = 5, dpi = 600)

cat("\n--- Caption text ---\n")
cat(sprintf(paste0("Leave-one-out analysis of the immunocompromised contrast. ",
                   "Each of the %d immunocompromised subjects was removed in ",
                   "turn and CellDMC was rerun, giving %d comparisons of %d ",
                   "versus %d subjects. Counts are plotted with a +1 offset on ",
                   "a log scale; dashed lines mark the full-cohort counts of ",
                   "%s. Removing a single subject increases the count by up to ",
                   "two orders of magnitude, and the subject responsible ",
                   "differs by compartment (%s), indicating that the observed ",
                   "counts depend on individual specimens rather than a stable ",
                   "cohort-level signal.\n"),
            length(ic_idx), length(ic_idx), length(ic_idx) - 1L,
            sum(ic_pheno_s3 == 0L),
            paste(sprintf("%s = %d", comp_labels[comps], obs_counts_s3[comps]),
                  collapse = ", "),
            driver_txt))


# -----------------------------------------------------------------------------
# V.7 Table 6: CellDMC counts and validation analyses
# -----------------------------------------------------------------------------
# One row per analysis, one column per compartment, so the observed counts and
# the three robustness checks read side by side.

ct_levels_t6  <- c("Epi", "Neutro", "Rest")
ct_display_t6 <- COMPARTMENT_LABELS
fmt_n <- function(x) format(round(as.numeric(x)), big.mark = ",", trim = TRUE)

loo_range_t6 <- vapply(ct_levels_t6, function(ct)
  sprintf("%s to %s", fmt_n(min(loo_tbl[[ct]], na.rm = TRUE)),
          fmt_n(max(loo_tbl[[ct]], na.rm = TRUE))), character(1))

perm_row <- function(col, fmt = fmt_n) {
  if (is.null(perm_summary_s3)) return(rep("\u2014", length(ct_levels_t6)))
  vapply(ct_levels_t6, function(ct) {
    v <- perm_summary_s3[[col]][perm_summary_s3$cell_type == ct]
    if (length(v) == 0) "\u2014" else fmt(v)
  }, character(1))
}

# The permutation null is the one object here that can arrive from a previous
# run, so it is the one that can silently disagree with the counts computed in
# THIS run -- which is exactly what happened before this guard existed: a null
# computed against observed counts of 118 / 235 / 83 was reported alongside
# freshly computed counts of 125 / 235 / 88, and the empirical p-values in
# Table 6 were therefore measured against the wrong anchors. The stopifnot in
# the cluster-load branch above covers only one of the ways perm_summary_s3 can
# be populated; this covers all of them, at the point of use.
if (!is.null(perm_summary_s3)) {
  .anchor <- obs_counts_s3[perm_summary_s3$cell_type]
  if (!isTRUE(all.equal(as.integer(perm_summary_s3$observed),
                        as.integer(.anchor)))) {
    stop("The permutation null was computed against different observed counts ",
         "than this run produced.\n",
         "  null was anchored to: ",
         paste(sprintf("%s = %d", perm_summary_s3$cell_type,
                       as.integer(perm_summary_s3$observed)), collapse = ", "),
         "\n  this run computed:    ",
         paste(sprintf("%s = %d", names(obs_counts_s3),
                       as.integer(obs_counts_s3)), collapse = ", "),
         "\n  The empirical p-values would be measured against the wrong ",
         "counts.\n  Recompute the null, or recompute emp_p_value from ",
         "perm_mat against the\n  current observed counts, before building ",
         "Table 6.")
  }
  rm(.anchor)
}

table6_celldmc <- tibble(
  Analysis = c("DMCs, unadjusted",
               "DMCs, adjusted",
               "DMCs, leave-one-out range",
               "Permuted median",
               "Permuted 95th percentile",
               "Empirical p"))
for (ct in ct_levels_t6) {
  i <- match(ct, ct_levels_t6)
  table6_celldmc[[unname(ct_display_t6[ct])]] <- c(
    fmt_n(obs_counts_s3[ct]),
    fmt_n(adj_counts_s3[ct]),
    unname(loo_range_t6[i]),
    perm_row("perm_median")[i],
    perm_row("perm_p95")[i],
    perm_row("emp_p_value", function(v) sprintf("%.3f", v))[i])
}

cat("\nTable 6:\n"); print(as.data.frame(table6_celldmc))

save_pub_table(
  table6_celldmc,
  file.path(TBL_DIR, "Table6_CellDMC_immunocompromised"),
  caption = paste0("Table 6. Cell-type-specific differentially methylated ",
                   "cytosines in the immunocompromised contrast and validation ",
                   "analyses."),
  note = paste0("Counts are CpGs reaching a Benjamini-Hochberg false discovery ",
                "rate below 0.05 in the contrast of ", sum(ic_pheno_s3 == 1L),
                " immunocompromised versus ", sum(ic_pheno_s3 == 0L),
                " non-immunocompromised subjects. The adjusted analysis added ",
                "age and array slide as model covariates. The leave-one-out ",
                "range spans the ", nrow(loo_tbl), " re-analyses in which each ",
                "immunocompromised subject was removed in turn. The ",
                "permutation null shuffled the immunocompromised label ",
                ifelse(is.null(perm_summary_s3), "B",
                       format(nrow(perm_ok_s3), big.mark = ",")),
                " times preserving group sizes; the empirical p-value is the ",
                "proportion of permutations whose count equaled or exceeded ",
                "the observed count, with add-one smoothing. Others is the ",
                "pooled non-epithelial, non-neutrophil compartment."))


# -----------------------------------------------------------------------------
# V.8 Checkpoint and text-ready summary
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
saveRDS(list(invariance_check = invariance_check,
             obs_counts   = obs_counts_s3,
             adj_counts   = adj_counts_s3,
             loo_tbl      = loo_tbl,
             perm_summary = perm_summary_s3,
             ic_pheno     = ic_pheno_s3,
             frac_cdmc    = frac_cdmc_s3),
        file.path(RDS_DIR, "checkpoint_celldmc_S3.rds"))

cat("\n--- Summary for the response letter (comment 1.6) ---\n")
cat(sprintf(
  paste0("The immunocompromised contrast compared %d immunocompromised with %d ",
         "non-immunocompromised subjects; PARDS vs control and survivor vs ",
         "non-survivor were not tested under the four-subject rule. The ",
         "near-invariance rule removed %d of the %d compartments. Unadjusted ",
         "counts were %s; with age and slide as covariates, %s.\n"),
  sum(ic_pheno_s3 == 1L), sum(ic_pheno_s3 == 0L),
  n_drop_methods, ncol(frac_cdmc_s3),
  paste(sprintf("%s = %d", names(obs_counts_s3), obs_counts_s3), collapse = ", "),
  paste(sprintf("%s = %d", names(adj_counts_s3), adj_counts_s3), collapse = ", ")))

if (!is.null(perm_summary_s3))
  cat(sprintf(
    paste0("Across %d label permutations preserving group sizes, median ",
           "permuted counts were %s against observed %s, giving empirical ",
           "p-values of %s.\n"),
    nrow(perm_ok_s3),
    paste(perm_summary_s3$perm_median, collapse = " / "),
    paste(perm_summary_s3$observed, collapse = " / "),
    paste(sprintf("%.4f", perm_summary_s3$emp_p_value), collapse = " / ")))

cat("\nPART V complete:", format(Sys.time()), "\n")




# =============================================================================
# PART VI. SENSITIVITY TO THE DECONVOLUTION AND ADJUSTMENT STRATEGY
# =============================================================================


# -----------------------------------------------------------------------------
# VI.1 Five adjustment strategies against no adjustment
# -----------------------------------------------------------------------------

# Five strategies spanning no adjustment to full adjustment. The prespecified
# strategy is the rule stated in the manuscript methods; the others are
# reported to show how the clustering result depends on that choice.

build_covars <- function(p, strategy) {
  ref <- cell_cols[which.max(colMeans(p))]
  switch(strategy,
    "Unadjusted" = NULL,
    "Prespecified (reference dropped + variance filter)" = {
      cand <- setdiff(cell_cols, ref)
      v    <- apply(p[, cand, drop = FALSE], 2, var)
      keep <- cand[v >= VAR_THRESHOLD]
      p[, keep, drop = FALSE]
    },
    "Epithelial and neutrophil only" = p[, c("Epi", "Neutro"), drop = FALSE],
    "All cell types except reference (no variance filter)" =
      p[, setdiff(cell_cols, ref), drop = FALSE],
    "Immune fraction only" = {
      m <- matrix(rowSums(p[, immune_cols, drop = FALSE]), ncol = 1,
                  dimnames = list(rownames(p), "Immune"))
      m
    })
}

# Computation order only. Display order for the reported table is set by
# STRATEGY_ORDER in S7.8, so the console print below is not in table order.
STRATEGIES <- c("Unadjusted",
                "Prespecified (reference dropped + variance filter)",
                "Epithelial and neutrophil only",
                "All cell types except reference (no variance filter)",
                "Immune fraction only")

adj_sensitivity <- map_dfr(ANALYTIC_SETS, function(s) {
  map_dfr(STRATEGIES, function(st) {
    cv <- build_covars(s$prop, st)
    mm <- if (is.null(cv)) s$m else removeBatchEffect(s$m, covariates = cv)
    ok <- rowSums(is.na(mm) | !is.finite(mm)) == 0
    mm <- mm[ok, , drop = FALSE]

    pv    <- rowVars(mm, na.rm = TRUE)
    mtop  <- mm[pv >= quantile(pv, 0.95, na.rm = TRUE), , drop = FALSE]
    pca   <- prcomp(t(mtop), center = TRUE, scale. = FALSE)
    ve    <- (pca$sdev^2) / sum(pca$sdev^2)
    d     <- dist(t(mtop)); hc <- hclust(d, method = "ward.D2")
    kmax  <- min(6, ncol(mm) - 1)
    sil   <- vapply(2:kmax, function(k) mean(silhouette(cutree(hc, k), d)[, 3]),
                    numeric(1))
    sizes2 <- table(cutree(hc, 2))

    tibble(analytic_set = s$label, strategy = st,
           n_covariates = if (is.null(cv)) 0L else ncol(cv),
           residual_df  = ncol(s$m) - (if (is.null(cv)) 0L else ncol(cv)) - 1L,
           probes_used  = nrow(mtop),
           pc1_pct      = ve[1] * 100,
           sil_k2       = sil[1],
           best_k       = (2:kmax)[which.max(sil)],
           best_sil     = max(sil),
           k2_sizes     = paste(sort(as.integer(sizes2), decreasing = TRUE),
                                collapse = " / "))
  })
})

cat("\nSensitivity to alternative adjustment strategies:\n")
print(as.data.frame(adj_sensitivity %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
write_csv(adj_sensitivity,
          file.path(SENS_TBL_DIR, "decon_adjustment_sensitivity.csv"))

# No publication table is written here. These five strategies are superseded by
# the six-strategy table assembled in S7.8, which adds the stage-one
# decomposition and is the version reported in the manuscript. Writing both
# left two near-identical documents in SENS_TBL_DIR, differing only by one row
# and by strategy order. The CSV above is retained as the intermediate.

# -----------------------------------------------------------------------------
# VI.2 How much of the methylome does composition explain (overadjustment)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# removeBatchEffect returns residuals, so the per-probe proportion of variance
# attributable to composition is 1 - var(adjusted) / var(unadjusted). This
# quantifies what the adjustment removes, which is the substance of the
# reviewer's overadjustment concern: in PARDS, epithelial loss and neutrophil
# influx may themselves be disease biology rather than technical noise.

overadj <- map_dfr(ANALYTIC_SETS, function(s) {
  cv <- build_covars(s$prop, "Prespecified (reference dropped + variance filter)")
  mm <- removeBatchEffect(s$m, covariates = cv)
  ok <- rowSums(is.na(mm) | !is.finite(mm)) == 0
  v_raw <- rowVars(s$m[ok, , drop = FALSE])
  v_adj <- rowVars(mm[ok, , drop = FALSE])
  r2    <- 1 - (v_adj / v_raw)
  r2    <- r2[is.finite(r2)]

  pca_raw <- prcomp(t(s$m[ok, , drop = FALSE][
    rowVars(s$m[ok, , drop = FALSE]) >=
      quantile(rowVars(s$m[ok, , drop = FALSE]), 0.95), , drop = FALSE]),
    center = TRUE, scale. = FALSE)
  pc1 <- pca_raw$x[, 1]

  tibble(analytic_set = s$label,
         median_var_explained = median(r2),
         q25 = quantile(r2, 0.25), q75 = quantile(r2, 0.75),
         pct_probes_over_50 = 100 * mean(r2 > 0.5),
         cor_Epi_PC1_raw    = cor(s$prop[names(pc1), "Epi"], pc1),
         cor_Neutro_PC1_raw = cor(s$prop[names(pc1), "Neutro"], pc1))
})

cat("\nVariance in the unadjusted methylome attributable to composition,\n")
cat("and correlation of the dominant fractions with unadjusted PC1:\n")
print(as.data.frame(overadj %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
write_csv(overadj, file.path(SENS_TBL_DIR, "decon_variance_explained.csv"))

cat("\nNOTE: the sign of the neutrophil-PC1 correlation printed above is the\n")
cat("value to quote in the Results (Reviewer 3, comment 12). PC axis sign is\n")
cat("arbitrary, so report it together with the epithelial correlation, which\n")
cat("has the opposite sign, rather than in isolation.\n")

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# VI.3 Stage-one (epithelial / fibroblast / immune) adjustment
# -----------------------------------------------------------------------------

# =============================================================================
# The seven-immune-subtype model rests on the 131 stage-two reference CpGs
# retained on EPIC v2, and S6 showed those subtype estimates fit poorly in
# epithelium-dominated specimens. The stage-one decomposition is better
# supported: more retained reference CpGs (Supplemental Table 4), three
# compartments, and adjustment on two
# covariates rather than seven, which matters at n = 13 to 24.
#
# This block re-derives the composition adjustment from stage one alone
# (epidish with centEpiFibIC), then repeats the variance filter, PCA, Ward.D2
# clustering, and silhouette in each of the three analytic sets.
#
# NOTE ON INTERPRETATION: hepidish partitions the stage-one immune fraction
# using stage two, so the immune total here equals the sum of the seven
# hepidish immune columns by construction. This analysis is therefore a
# better-justified version of the "immune fraction only" strategy in S6.4, not
# an independent estimate. The block verifies that identity explicitly.

cat("\n=== PART VI.3. Stage-one (Epi / Fib / IC) adjustment sensitivity ===\n")

immune_cols_s7 <- setdiff(cell_cols, c("Epi", "Fib"))

# Partitions from the prespecified analyses, for concordance comparison
prespec_partition <- list(
  "Full cohort (n = 24)" = setNames(sample_sheet$Adjusted_Group,
                                    sample_sheet$Sample_Name),
  "PARDS only"           = setNames(assign_pards$PARDS_only_Group,
                                    assign_pards$Sample_Name),
  "New PARDS only"       = setNames(paste0("NewGroup", cl_new), names(cl_new))
)

run_stage1 <- function(s, prespec_key) {
  cat("\n---", s$label, "---\n")

  # Stage-one estimates were computed once in section III.3 and are reused
  # here so the log ratio in Table 3 and the adjustment below cannot diverge.
  e1 <- stage1_by_set[[s$label]]
  cat("Stage-one reference CpGs used:",
      sum(rownames(s$beta) %in% rownames(centEpiFibIC.m)), "\n")

  # Verify the identity against the hepidish output
  ic_hep <- rowSums(s$prop[, immune_cols_s7, drop = FALSE])
  cat("Max |IC(stage 1) - sum(hepidish immune)|:",
      signif(max(abs(e1[, "IC"] - ic_hep[rownames(e1)])), 3), "\n")

  comp_pct <- colMeans(e1) * 100
  cat("Mean compartment fractions (%): ",
      paste(sprintf("%s = %.1f", names(comp_pct), comp_pct), collapse = ", "),
      "\n", sep = "")

  # Same two prespecified reductions, applied to the three compartments
  ref1 <- colnames(e1)[which.max(colMeans(e1))]
  cand <- setdiff(colnames(e1), ref1)
  v    <- apply(e1[, cand, drop = FALSE], 2, var)
  keep <- cand[v >= VAR_THRESHOLD]
  cat("Reference compartment dropped:", ref1,
      "| covariates retained:", paste(keep, collapse = ", "),
      "| dropped for low variance:",
      ifelse(length(setdiff(cand, keep)) == 0, "none",
             paste(setdiff(cand, keep), collapse = ", ")), "\n")

  cov1 <- e1[, keep, drop = FALSE]
  mm   <- removeBatchEffect(s$m, covariates = cov1)
  ok   <- rowSums(is.na(mm) | !is.finite(mm)) == 0
  mm   <- mm[ok, , drop = FALSE]
  rdf  <- ncol(s$m) - ncol(cov1) - 1
  cat("Residual degrees of freedom:", rdf, "\n")

  res <- cluster_diagnostics(mm, paste0(s$label, " [stage-one adjustment]"))

  k2   <- min(2, ncol(mm) - 1)
  cl2  <- cutree(res$hc, k = 2)
  sizes <- paste(sort(as.integer(table(cl2)), decreasing = TRUE), collapse = " / ")
  cat("Cluster sizes at k = 2:", sizes, "\n")

  # Concordance with the prespecified seven-covariate partition
  pp  <- prespec_partition[[prespec_key]]
  common <- intersect(names(cl2), names(pp))
  ari <- adjusted_rand_index(pp[common], paste0("S1Group", cl2[common]))
  cat("Adjusted Rand index vs the prespecified partition:", round(ari, 3), "\n")

  # Variance in the unadjusted methylome removed by this adjustment
  r2 <- 1 - (rowVars(mm) / rowVars(s$m[ok, , drop = FALSE]))
  r2 <- r2[is.finite(r2)]

  list(
    summary = tibble(
      analytic_set = s$label,
      n_covariates = ncol(cov1),
      covariates   = paste(keep, collapse = " + "),
      residual_df  = rdf,
      pc1_pct      = res$var_exp[1] * 100,
      sil_k2       = res$sil$mean_silhouette[res$sil$k == 2],
      best_k       = res$best_k,
      best_sil     = res$best_sil,
      k2_sizes     = sizes,
      ari_vs_prespecified = ari,
      median_var_removed  = median(r2)),
    sil = res$sil %>% mutate(analytic_set = s$label),
    est = as.data.frame(e1) %>% rownames_to_column("Sample_Name") %>%
      mutate(analytic_set = s$label),
    res = res, cl2 = cl2)
}

s7_full <- run_stage1(ANALYTIC_SETS[[1]], "Full cohort (n = 24)")
s7_pards <- run_stage1(ANALYTIC_SETS[[2]], "PARDS only")
s7_new  <- run_stage1(ANALYTIC_SETS[[3]], "New PARDS only")

s7_summary <- bind_rows(s7_full$summary, s7_pards$summary, s7_new$summary)
cat("\nStage-one adjustment summary:\n")
print(as.data.frame(s7_summary %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
write_csv(s7_summary, file.path(SENS_TBL_DIR, "stage1_adjustment_summary.csv"))

s7_sil <- bind_rows(s7_full$sil, s7_pards$sil, s7_new$sil)
write_csv(s7_sil, file.path(SENS_TBL_DIR, "stage1_adjustment_silhouette_by_k.csv"))

s7_est <- bind_rows(s7_full$est, s7_pards$est, s7_new$est) %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID), by = "Sample_Name") %>%
  dplyr::select(analytic_set, Subject_ID, Epi, Fib, IC)
write_csv(s7_est, file.path(SENS_TBL_DIR, "stage1_compartment_estimates.csv"))

save_support_table(
  s7_summary %>%
    transmute(`Analytic set` = analytic_set,
              `Covariates` = covariates,
              `Covariates (n)` = n_covariates,
              `Residual df` = residual_df,
              `PC1 variance (%)` = sprintf("%.1f", pc1_pct),
              `Silhouette at k = 2` = sprintf("%.3f", sil_k2),
              `Cluster sizes at k = 2` = k2_sizes,
              `ARI vs prespecified` = sprintf("%.3f", ari_vs_prespecified)),
  file.path(SENS_TBL_DIR, "Supporting_stage1_adjustment"),
  caption = paste0("Supporting table. Composition adjustment using the ",
                   "stage-one epithelial, fibroblast, and immune cell ",
                   "decomposition only."),
  note = paste0("The stage-one decomposition uses the centEpiFibIC reference, ",
                "of which ", ref_overlap$ref1_retained[1], " of ",
                ref_overlap$ref1_total[1], " CpGs were retained on the EPIC v2 ",
                "array, compared with ", ref_overlap$ref2_retained[1], " of ",
                ref_overlap$ref2_total[1], " for the stage-two immune subtype ",
                "reference. Adjusting on the stage-one compartments requires ",
                "fewer covariates and is therefore less prone to overfitting ",
                "at these sample sizes. The same two prespecified reductions ",
                "were applied: the compartment with the largest mean fraction ",
                "was dropped as the reference, and any remaining compartment ",
                "with variance below 1e-4 was dropped. The adjusted Rand index ",
                "compares each partition with that obtained under the ",
                "prespecified seven-covariate adjustment. Cluster sizes are ",
                "reported alongside silhouette width because a cluster of one ",
                "specimen inflates the mean without indicating structure."))

# (No figure is written; Table 7 reports the stage-one adjustment beside
#  the other five strategies.)

saveRDS(list(summary = s7_summary, sil = s7_sil, estimates = s7_est,
             full = s7_full, pards = s7_pards, new = s7_new),
        file.path(RDS_DIR, "checkpoint_stage1_adjustment_S7.rds"))



# -----------------------------------------------------------------------------
# VI.4 Table 7: the six strategies together
# -----------------------------------------------------------------------------

# S6.4 produced five strategies and S7 added the stage-one decomposition. The
# response letter and the Results describe them as one comparison, so they are
# combined into a single supplementary table here, ordered to follow the
# deconvolution hierarchy.
#
# This table supersedes SuppTable_adjustment_sensitivity written in S6.4, which
# is no longer rendered.

# Ordered as a progression through the deconvolution hierarchy rather than by
# model size: no adjustment, the stage-one decomposition, its immune
# compartment alone, the prespecified rule combining both stages, the same rule
# without the variance filter, and finally the two dominant cell types alone.
# Note that the covariate count is therefore not monotonic down the column
# (0, 2, 1, 7, 8, 2), which is intended.
STRATEGY_ORDER <- c(
  "Unadjusted",
  "Stage-one decomposition (epithelial, fibroblast, immune)",
  "Immune fraction only",
  "Prespecified (reference dropped + variance filter)",
  "All cell types except reference (no variance filter)",
  "Epithelial and neutrophil only")

s7_rows <- s7_summary %>%
  transmute(analytic_set,
            strategy = "Stage-one decomposition (epithelial, fibroblast, immune)",
            n_covariates, residual_df, pc1_pct, sil_k2, best_k, best_sil,
            k2_sizes)

adj_all <- adj_sensitivity %>%
  dplyr::select(analytic_set, strategy, n_covariates, residual_df,
                pc1_pct, sil_k2, best_k, best_sil, k2_sizes) %>%
  bind_rows(s7_rows) %>%
  mutate(strategy = factor(strategy, levels = STRATEGY_ORDER),
         analytic_set = factor(analytic_set,
                               levels = vapply(ANALYTIC_SETS, `[[`,
                                               character(1), "label"))) %>%
  arrange(analytic_set, strategy)

cat("\nMerged adjustment sensitivity (six strategies):\n")
print(as.data.frame(adj_all %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
write_csv(adj_all,
          file.path(SENS_TBL_DIR, "decon_adjustment_sensitivity_all.csv"))


# Abbreviated labels, written so the progression through the deconvolution
# hierarchy is visible in the table itself. Lookup is by name, so the order
# here is for readability only; display order is set by STRATEGY_ORDER above.
STRATEGY_SHORT <- c(
  "Unadjusted"                                               = "Unadjusted",
  "Stage-one decomposition (epithelial, fibroblast, immune)" = "Stage one: epithelial, fibroblast, immune",
  "Immune fraction only"                                     = "Stage one: immune fraction only",
  "Prespecified (reference dropped + variance filter)"        = "Prespecified: stages one and two, variance filter",
  "All cell types except reference (no variance filter)"     = "Stages one and two, no variance filter",
  "Epithelial and neutrophil only"                           = "Epithelial and neutrophil only")

adj_all <- adj_all %>%
  mutate(strategy = factor(as.character(strategy), levels = STRATEGY_ORDER)) %>%
  arrange(analytic_set, strategy)
stopifnot(!any(is.na(adj_all$strategy)))

table7_sensitivity <- adj_all %>%
  transmute(Cohort                 = as.character(analytic_set),
            `Adjustment strategy`  = unname(STRATEGY_SHORT[as.character(strategy)]),
            `Covariates (n)`       = as.character(n_covariates),
            `Residual df`          = as.character(residual_df),
            `PC1 (%)`              = sprintf("%.1f", pc1_pct),
            `Silhouette, k = 2`    = sprintf("%.3f", sil_k2),
            `Cluster sizes, k = 2` = k2_sizes)
stopifnot(!any(is.na(table7_sensitivity$`Adjustment strategy`)))

save_grouped_table(
  table7_sensitivity, "Cohort",
  file.path(TBL_DIR, "Table7_adjustment_sensitivity"),
  caption = paste0("Table 7. Sensitivity of the composition adjustment and ",
                   "resulting cluster structure to six adjustment strategies."),
  note = paste0(
    "For each analytic set and strategy, cell-type proportions were supplied ",
    "to limma::removeBatchEffect as continuous covariates, the top 5% most ",
    "variable CpGs were reselected, and Ward.D2 hierarchical clustering of ",
    "euclidean distances was recomputed. Strategies are ordered to follow the ",
    "deconvolution hierarchy: no adjustment, the stage-one decomposition, its ",
    "immune compartment alone, the prespecified rule combining stages one and ",
    "two, the same rule without the variance filter, and adjustment for the ",
    "two dominant cell types alone. Strategies are abbreviated as follows. ",
    "Stage one, the epithelial, fibroblast, and immune compartments from the ",
    "centEpiFibIC reference. Stage one immune fraction only, the summed immune ",
    "compartment as a single covariate. Prespecified, the cell type with the ",
    "largest mean fraction dropped as the reference and any remaining fraction ",
    "with variance below 1e-4 dropped. Epithelial and neutrophil only, those ",
    "two fractions, which are the only strategy not derived from the ",
    "deconvolution hierarchy. Of the reference CpGs, ",
    ref_overlap$ref1_retained[1], " of ", ref_overlap$ref1_total[1],
    " were retained on the EPIC v2 array for the stage-one reference and ",
    ref_overlap$ref2_retained[1], " of ", ref_overlap$ref2_total[1],
    " for the stage-two immune subtype reference; because the three stage-one ",
    "compartments sum to one, dropping any one as the reference yields ",
    "identical residuals. In the ", n_new, "-specimen set the prespecified and ",
    "no-variance-filter strategies give identical results because the only ",
    "covariate removed by the variance filter is the eosinophil fraction, ",
    "which is estimated at zero in nearly all specimens and therefore ",
    "contributes no adjustment. Cluster sizes are reported alongside the ",
    "silhouette width because a cluster containing one or two specimens ",
    "inflates the mean silhouette without indicating subgroup structure."),
  widths = c(2.35, 0.80, 0.70, 0.70, 0.90, 0.85))

saveRDS(adj_all, file.path(RDS_DIR, "checkpoint_adjustment_sensitivity_all.rds"))
cat("\nPART VI complete:", format(Sys.time()), "\n")



# =============================================================================
# PART VII. RESPONSE-LETTER ANALYSES AND SESSION INFORMATION
# =============================================================================

# The analysis below is reported in the point-by-point response to Reviewer 1,
# comment 5, and is not part of the manuscript Results. It is retained here so
# the repository reproduces every number given to the reviewers.

# =============================================================================
# Reviewer 1 asks whether cell-type deconvolution invalidates the earlier work
# (Williams et al., Respir Res 2022;23:181). Rather than reprocessing that
# dataset, we ask a narrower and directly answerable question: does the nasal
# cell composition estimated in the current cohort explain the nasal
# transcriptomic subgroup assignment (Subgroups A to D) reported in 2022?
#
# Design decisions, made to avoid repeating the problems reviewers flagged
# elsewhere in this manuscript:
#   - ONE omnibus multivariate test (PERMANOVA) on the full composition vector,
#     rather than a univariate test per cell type. Subgroups B and D contain
#     two subjects each, so four-group univariate contrasts would reproduce
#     exactly the small-n comparisons objected to in comments 1.6 and 2.4.
#   - R-squared is the quantity of interest (how much compositional variance
#     subgroup explains), not the p-value alone. With 16 subjects and 3 model
#     degrees of freedom, R-squared is upward biased; the permutation p carries
#     the inference.
#   - betadisper() is run alongside, because PERMANOVA can be driven by unequal
#     within-group dispersion rather than by differences in group centroids.
#     Within-subgroup spread is very large (Subgroup C spans 0.2% to 100%
#     epithelial fraction), so this check is not optional.
#   - Bray-Curtis on the raw proportions is used as the primary distance. It
#     tolerates the exact zeros in these estimates, unlike a CLR or Aitchison
#     approach, which would require an arbitrary zero-replacement.
#
# Transcriptomic subgroup assignments were read from Fig. 4A of the 2022 paper,
# taking the day 1 assignment where available and the earliest available
# specimen otherwise. ARDS_003 (day 7) and ARDS_023 (day 14) are the only two
# not from day 1; a sensitivity analysis excludes them. The specimens are the
# same brushings analyzed in both studies; the methylomic and transcriptomic
# assays were simply run on different schedules.

cat("\n=== PART VII.1. Composition vs 2022 transcriptomic subgroup ===\n")

transcriptomic_subgroup <- tibble::tribble(
  ~Subject_ID,  ~Transcriptomic_Subgroup, ~Assignment_Day,
  "ARDS_001",   "C",                      "Day 1",
  "ARDS_002",   "C",                      "Day 1",
  "ARDS_003",   "D",                      "Day 7",
  "ARDS_004",   "B",                      "Day 1",
  "ARDS_010",   "C",                      "Day 1",
  "ARDS_015",   "C",                      "Day 1",
  "ARDS_016",   "A",                      "Day 1",
  "ARDS_018",   "A",                      "Day 1",
  "ARDS_019",   "A",                      "Day 1",
  "ARDS_021",   "C",                      "Day 1",
  "ARDS_023",   "A",                      "Day 14",
  "ARDS_026",   "A",                      "Day 1",
  "ARDS_029",   "D",                      "Day 1",
  "ARDS_030",   "B",                      "Day 1",
  "ARDS_031",   "C",                      "Day 1",
  "ARDS_032",   "A",                      "Day 1"
)

write_csv(transcriptomic_subgroup,
          file.path(SENS_TBL_DIR, "transcriptomic_subgroup_assignments_2022.csv"))

# -----------------------------------------------------------------------------
# S2.1 Assemble the analytic set
# -----------------------------------------------------------------------------
comp_tx <- cell_props_pards %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID), by = "Sample_Name") %>%
  inner_join(transcriptomic_subgroup, by = "Subject_ID") %>%
  mutate(
    Subgroup    = factor(Transcriptomic_Subgroup, levels = c("A", "B", "C", "D")),
    Subgroup_BD = factor(if_else(Transcriptomic_Subgroup %in% c("B", "D"),
                                 "B or D", "A or C"))
  ) %>%
  arrange(Subgroup, Subject_ID)

cat("Subjects with both composition estimates and a 2022 subgroup:",
    nrow(comp_tx), "of", n_pards, "\n")
print(table(comp_tx$Subgroup, dnn = "Transcriptomic subgroup"))
cat("Assignments not taken from day 1:",
    paste(comp_tx$Subject_ID[comp_tx$Assignment_Day != "Day 1"], collapse = ", "), "\n")

comp_mat <- as.matrix(comp_tx[, cell_cols, drop = FALSE])
rownames(comp_mat) <- comp_tx$Subject_ID

# -----------------------------------------------------------------------------
# S2.2 Descriptive composition by subgroup (no univariate p-values by design)
# -----------------------------------------------------------------------------
other_immune <- intersect(c("B", "NK", "CD4T", "CD8T", "Mono", "Eosino"), cell_cols)

fmt_med_iqr <- function(x) sprintf("%.1f [%.1f, %.1f]",
                                   median(x) * 100,
                                   quantile(x, 0.25) * 100,
                                   quantile(x, 0.75) * 100)

desc_by_subgroup <- comp_tx %>%
  mutate(Other_immune = rowSums(across(all_of(other_immune)))) %>%
  group_by(Subgroup) %>%
  summarise(
    n              = n(),
    Epithelial     = fmt_med_iqr(Epi),
    Fibroblast     = fmt_med_iqr(Fib),
    Neutrophil     = fmt_med_iqr(Neutro),
    `Other immune` = fmt_med_iqr(Other_immune),
    .groups = "drop"
  )

cat("\nCell composition by transcriptomic subgroup, median [IQR] percent:\n")
print(as.data.frame(desc_by_subgroup))

save_support_table(
  desc_by_subgroup %>% rename(`Transcriptomic subgroup` = Subgroup),
  file.path(SENS_TBL_DIR, "Supporting_composition_by_transcriptomic_subgroup"),
  caption = paste0("Supporting table. Estimated nasal cell composition by ",
                   "nasal transcriptomic subgroup previously reported in ",
                   "Williams et al., Respir Res 2022;23:181."),
  note = paste0("Values are median [interquartile range] percent of the ",
                "EpiDISH (hepidish, RPC) composition estimate for the ",
                nrow(comp_tx), " subjects with both a current methylation ",
                "specimen and a transcriptomic subgroup assignment. Other ",
                "immune sums B, NK, CD4+ T, CD8+ T, monocyte, and eosinophil ",
                "fractions. Subgroups B and D contain two subjects each; ",
                "univariate between-subgroup p-values are deliberately not ",
                "reported. Association was assessed by a single multivariate ",
                "permutational test (see PERMANOVA supplementary table)."))

# -----------------------------------------------------------------------------
# S2.3 PERMANOVA: does subgroup explain compositional variance?
# -----------------------------------------------------------------------------
set.seed(20260101)

d_bray <- vegdist(comp_mat, method = "bray")

pmv_4 <- adonis2(d_bray ~ Subgroup, data = comp_tx, permutations = N_PERM)
cat("\nPERMANOVA, Bray-Curtis, four subgroups (n =", nrow(comp_tx), "):\n")
print(pmv_4)

pmv_bd <- adonis2(d_bray ~ Subgroup_BD, data = comp_tx, permutations = N_PERM)
cat("\nPERMANOVA, Bray-Curtis, B or D vs A or C (the 2022 primary contrast):\n")
print(pmv_bd)

# Sensitivity 1: Euclidean distance instead of Bray-Curtis
d_euc   <- dist(comp_mat, method = "euclidean")
pmv_euc <- adonis2(d_euc ~ Subgroup, data = comp_tx, permutations = N_PERM)
cat("\nPERMANOVA, Euclidean distance, four subgroups:\n")
print(pmv_euc)

# Sensitivity 2: day 1 assignments only
comp_d1 <- comp_tx %>% filter(Assignment_Day == "Day 1")
mat_d1  <- as.matrix(comp_d1[, cell_cols, drop = FALSE])
rownames(mat_d1) <- comp_d1$Subject_ID
pmv_d1 <- adonis2(vegdist(mat_d1, method = "bray") ~ Subgroup,
                  data = comp_d1, permutations = N_PERM)
cat("\nPERMANOVA, Bray-Curtis, day 1 assignments only (n =", nrow(comp_d1), "):\n")
print(pmv_d1)

# -----------------------------------------------------------------------------
# S2.4 Dispersion check
# -----------------------------------------------------------------------------
bd_disp   <- betadisper(d_bray, comp_tx$Subgroup)
disp_test <- permutest(bd_disp, permutations = N_PERM)
cat("\nHomogeneity of multivariate dispersion (betadisper):\n")
print(disp_test)
cat("\nMean distance to subgroup centroid:\n")
print(round(tapply(bd_disp$distances, comp_tx$Subgroup, mean), 3))

# -----------------------------------------------------------------------------
# S2.5 Results table
# -----------------------------------------------------------------------------
pull_pmv <- function(p, label) {
  tibble(Analysis   = label,
         `R2`       = sprintf("%.3f", p$R2[1]),
         `F`        = sprintf("%.2f", p$F[1]),
         `p (perm)` = sprintf("%.3f", p$`Pr(>F)`[1]))
}

permanova_table <- bind_rows(
  pull_pmv(pmv_4,   paste0("Four subgroups, Bray-Curtis (n = ", nrow(comp_tx), ")")),
  pull_pmv(pmv_bd,  "B or D vs A or C, Bray-Curtis"),
  pull_pmv(pmv_euc, "Four subgroups, Euclidean"),
  pull_pmv(pmv_d1,  paste0("Four subgroups, day 1 only (n = ", nrow(comp_d1), ")"))
)

cat("\nPERMANOVA summary:\n")
print(as.data.frame(permanova_table))
write_csv(permanova_table,
          file.path(SENS_TBL_DIR, "permanova_composition_vs_transcriptomic.csv"))

save_pub_table(
  permanova_table,
  file.path(SENS_TBL_DIR, "Supporting_PERMANOVA_composition_vs_transcriptomic"),
  caption = paste0("Supporting table. Permutational multivariate analysis of ",
                   "variance testing whether nasal transcriptomic subgroup ",
                   "explains estimated cell composition."),
  note = paste0("Each row is a single omnibus test across all nine estimated ",
                "cell-type fractions, with ", format(N_PERM, big.mark = ","),
                " permutations. R-squared is the proportion of compositional ",
                "variance attributable to subgroup. Within-group dispersion ",
                "was assessed separately by betadisper (permutation p = ",
                sprintf("%.3f", disp_test$tab$`Pr(>F)`[1]),
                "); a difference in dispersion rather than in group centroids ",
                "can also produce a low PERMANOVA p-value."))

# -----------------------------------------------------------------------------
# VII.1e Ordination of composition by transcriptomic subgroup
# -----------------------------------------------------------------------------
pcoa     <- cmdscale(d_bray, k = 2, eig = TRUE)
pcoa_var <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])
# Ordination retained for the checkpoint below; no figure is written, since
# the response to Reviewer 1 comment 5 quotes the PERMANOVA result only.

# -----------------------------------------------------------------------------
# S2.7 Checkpoint and text-ready summary
# -----------------------------------------------------------------------------
saveRDS(list(
  transcriptomic_subgroup = transcriptomic_subgroup, comp_tx = comp_tx,
  desc_by_subgroup = desc_by_subgroup, permanova_table = permanova_table,
  pmv_4 = pmv_4, pmv_bd = pmv_bd, pmv_euc = pmv_euc, pmv_d1 = pmv_d1,
  disp_test = disp_test, pcoa = pcoa
), file.path(RDS_DIR, "checkpoint_sensitivity_transcriptomic.rds"))

cat("\n--- Summary for the response letter (comment 1.5) ---\n")

# The dispersion clause is written from the betadisper result rather than
# hardcoded, so the sentence cannot contradict the test it cites. Note that a
# non-significant betadisper does not establish equal dispersion; with two
# subjects in two of the subgroups the test has very little power, so the
# wording claims only that no difference was detected.
disp_p      <- disp_test$tab$`Pr(>F)`[1]
disp_clause <- if (disp_p < 0.05) {
  sprintf(paste0("Within-subgroup dispersion differed (betadisper p = %.3f), so ",
                 "the PERMANOVA result should not be read as a comparison of ",
                 "group centroids alone"), disp_p)
} else {
  sprintf(paste0("No difference in within-subgroup dispersion was detected ",
                 "(betadisper p = %.3f), supporting interpretation of the ",
                 "PERMANOVA as a comparison of group centroids"), disp_p)
}

cat(sprintf(
  paste0("Among the %d subjects with both a current methylation specimen and a ",
         "nasal transcriptomic subgroup assignment from the 2022 report, ",
         "transcriptomic subgroup explained %.1f%% of the variance in estimated ",
         "nasal cell composition (PERMANOVA on Bray-Curtis distances, ",
         "%s permutations, F = %.2f, p = %.3f). %s. Restricting to day 1 ",
         "subgroup assignments (n = %d) gave R2 = %.3f, p = %.3f.\n"),
  nrow(comp_tx), pmv_4$R2[1] * 100, format(N_PERM, big.mark = ","),
  pmv_4$F[1], pmv_4$`Pr(>F)`[1], disp_clause,
  nrow(comp_d1), pmv_d1$R2[1], pmv_d1$`Pr(>F)`[1]))

cat("\nPART VII.1 complete:", format(Sys.time()), "\n")



# -----------------------------------------------------------------------------
# VII.2 Session information
# -----------------------------------------------------------------------------


cat("\nAnalysis complete:", format(Sys.time()), "\n")
cat("\nFigures written to: ", FIG_DIR, "\n", sep = "")
cat("Tables written to:  ", TBL_DIR, "\n", sep = "")
cat("Supporting output:  ", SENS_FIG_DIR, " and ", SENS_TBL_DIR, "\n", sep = "")

print(sessionInfo())

# Write the version manifest to its own file as well as to the log, so the exact
# package versions behind a given run can be cited or diffed without digging
# through a timestamped log.
# LOG_DIR, not PROJECT_ROOT: it is the one directory the script has already
# created and written to, so this cannot fail with the log sink still open.
si_path <- file.path(LOG_DIR, "sessionInfo.txt")
writeLines(c(paste("Run completed:", format(Sys.time())),
             "Script:       Methyl_PARDS_pipeline.R",
             "GEO series:   GSE337899",
             "",
             capture.output(sessionInfo())), si_path)
cat("Session information written to: ", si_path, "\n", sep = "")

sink()

