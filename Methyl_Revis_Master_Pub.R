# =============================================================================
# Pipeline overview:
#   1.   Setup and packages
#   2.   Load sample sheet + IDATs (sesame frontend for EPIC v2.0)
#   3-4. Per-sample QC: detection p-values, bisulfite conversion, intensity
#   5.   Beta and M-value matrices
#   6.   Probe filtering (cross-reactive, SNP-adjacent, detP, invariant)
#   7.   Collapse EPIC v2.0 replicate probes to base IDs
#   8.   Cell-type deconvolution (EpiDISH)
#   9.   Composition x raw-PCA / clinical associations
#   10.  Composition-adjusted M-values (removeBatchEffect)
#   11.  PCA on the adjusted matrix
#   12.  Hierarchical clustering with silhouette-based k selection
#   13.  Clinical-outcome analysis by adjusted groups
#   14.  Genomic annotation + enrichment setup
#   15.  Adjusted-group DMCs (limma) + promoter enrichment; CellDMC + validation
#   16.  Outputs and session info
#
#
# Requires R 4.5.0. Set PROJECT_ROOT (below) to your environment before running.

# =============================================================================


# -----------------------------------------------------------------------------
# 1. Setup, logging, and packages
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(sesame)
  library(sesameData)
  library(BiocParallel)  
  library(minfi)
  library(IlluminaHumanMethylationEPICv2manifest)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  library(limma)
  library(matrixStats)
  library(EpiDISH)
  library(pheatmap)
  library(RColorBrewer)
  library(cluster)
  library(clusterProfiler)
  library(msigdbr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(scatterplot3d)
  library(ggrepel)      
  library(patchwork)    
  library(wateRmelon)
})

options(ExperimentHub.ASK = FALSE)

# -----------------------------------------------------------------------------
# Route sesameData lookups to the alternate (Zhou-hosted) mirror
# -----------------------------------------------------------------------------
# Several EPICv2 sesameData objects (notably the KYCG.EPICv2 mask that pOOBAH
# pulls internally) are archived on ExperimentHub and fail to fetch. Setting
# SESAMEDATA_USE_ALT = TRUE resolves them from sesame's alternate mirror and
# caches them in-session. This must stay TRUE for the whole run.
options(SESAMEDATA_USE_ALT = TRUE)

set.seed(20260101)

# -----------------------------------------------------------------------------
# Configuration (edit PROJECT_ROOT to your environment)
# -----------------------------------------------------------------------------
PROJECT_ROOT <- "."   # <-- set to your project directory
IDAT_DIR     <- file.path(PROJECT_ROOT, "Combined IDATs")

CLINICAL     <- NULL  # clinical data already merged into the unified sheet
FIG_DIR      <- file.path(PROJECT_ROOT, "figures")
TBL_DIR      <- file.path(PROJECT_ROOT, "tables")
LOG_DIR      <- file.path(PROJECT_ROOT, "logs")
RDS_DIR      <- file.path(PROJECT_ROOT, "analysis_output")
for (d in c(FIG_DIR, TBL_DIR, LOG_DIR, RDS_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Publication-ready table writer
# -----------------------------------------------------------------------------

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

# Per-sample array QC thresholds (tune to your cohort before reporting)
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
# EPICv2 address file: local build for IDAT decoding
# -----------------------------------------------------------------------------
cat("Downloading EPICv2 manifest from the Zhou InfiniumAnnotation host...\n")
epicv2_tsv  <- sesameAnno_download("EPICv2.hg38.manifest.tsv.gz")
epicv2_addr <- sesameAnno_buildAddressFile(epicv2_tsv)
cat("  Address file built from EPICv2.hg38.manifest.tsv.gz\n")


# -----------------------------------------------------------------------------
# 2. Load sample sheet and iDATs
# -----------------------------------------------------------------------------


sample_sheet <- read_csv(SAMPLE_SHEET, show_col_types = FALSE) %>%
  mutate(Sex = na_if(Sex, "NA"))

cat("Loaded unified sample sheet:", nrow(sample_sheet), "samples,",
    ncol(sample_sheet), "columns\n")
cat("Columns:", paste(names(sample_sheet), collapse = ", "), "\n")
cat("Sample_Group counts:\n")
print(table(sample_sheet$Sample_Group, useNA = "ifany"))


# -----------------------------------------------------------------------------
# 3. Process iDATs with sesame openSesame (Noob + dye-bias)
# -----------------------------------------------------------------------------

idat_prefixes <- file.path(IDAT_DIR, sample_sheet$Sample_Name)

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
# 4. Per-sample detection p-values + QC
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
pdf(file.path(FIG_DIR, "qc_detP_per_sample.pdf"), width = 7, height = 4)
par(mar = c(8, 4, 2, 2))
barplot(mean_detP, las = 2,
        col = ifelse(mean_detP < 0.01, "steelblue", "firebrick"),
        ylab = "Mean detection P", main = "pOOBAH detection P by sample")
abline(h = 0.01, lty = 2)
dev.off()


# -----------------------------------------------------------------------------
# 4b. Extended per-sample QC metrics + supplementary figures
# -----------------------------------------------------------------------------

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
qc_extra <- map_dfr(sample_sheet$Sample_Name, function(s) {
  sdf <- sdf_list[[s]]
  gct <- tryCatch(bisConversionControl(sdf), error = function(e) NA_real_)
  mu  <- tryCatch(signalMU(sdf),            error = function(e) NULL)
  tibble(
    Sample_Name    = s,
    gct_score      = gct,
    median_meth    = if (!is.null(mu)) median(mu$M, na.rm = TRUE) else NA_real_,
    median_unmeth  = if (!is.null(mu)) median(mu$U, na.rm = TRUE) else NA_real_
  )
})

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
  
  # Heatmap of the concordance matrix, annotated by subject
  subj_ann <- data.frame(Subject = factor(subj_of[samp]))
  rownames(subj_ann) <- samp
  pheatmap(concord,
           annotation_row = subj_ann, annotation_col = subj_ann,
           display_numbers = (length(samp) <= 20), number_format = "%.2f",
           color = colorRampPalette(brewer.pal(9, "Blues"))(100),
           main = "EPIC rs-probe genotype concordance",
           filename = file.path(FIG_DIR, "SuppFig_genotype_identity.pdf"),
           width = 9, height = 8)
  pheatmap(concord,
           annotation_row = subj_ann, annotation_col = subj_ann,
           display_numbers = (length(samp) <= 20), number_format = "%.2f",
           color = colorRampPalette(brewer.pal(9, "Blues"))(100),
           main = "EPIC rs-probe genotype concordance",
           filename = file.path(FIG_DIR, "SuppFig_genotype_identity.png"),
           width = 9, height = 8)
  cat("  Identity outputs: genotype_identity_pairs.csv, SuppFig_genotype_identity.pdf/.png\n")
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
has_gct   <- any(!is.na(qc_plot_df$gct_score))
has_bscon <- any(!is.na(qc_plot_df$bs_conversion))
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
pInt <- ggplot(qc_plot_df, aes(log2_unmeth, log2_meth, color = status)) +
  # Actual exclusion rule: mean(log2 M, log2 U) < LOG2INT_MIN, which is the line
  # log2_meth = 2*LOG2INT_MIN - log2_unmeth. Samples below/left of it are flagged.
  geom_abline(slope = -1, intercept = 2 * LOG2INT_MIN,
              linetype = "dashed", color = "grey30") +
  geom_point(size = 2.4) +
  geom_text_repel(aes(label = label), size = 2.6, show.legend = FALSE,
                  min.segment.length = 0,        # draw a stem for EVERY label
                  segment.size = 0.3, segment.color = "grey55",
                  box.padding = 0.5, point.padding = 0.25,
                  max.overlaps = Inf, seed = 1) +
  scale_color_manual(values = status_cols, name = NULL) +
  labs(title = "Signal Intensity",
       x = expression(log[2]~median~unmethylated),
       y = expression(log[2]~median~methylated)) +
  theme_qc

# Sex concordance is intentionally not shown: three subjects have no recorded
# sex, so a predicted-vs-recorded concordance check is uninformative. The
# inferSex() prediction is retained in the QC table for covariate imputation.

supp_qc <- (pDetP | pBis) / pInt +
  plot_annotation(
    tag_levels = "A",
    title = "Supplementary Figure 1. Per-Sample Array Quality Control (EPIC v2.0)",
    subtitle = paste0("Flag thresholds: mean detP \u2264 ", DETP_SAMPLE_MAX,
                      ", GCT \u2264 ", GCT_MAX,
                      ", bisulfite conversion \u2265 ", BSCON_MIN, "%",
                      ", mean log2 intensity \u2265 ", LOG2INT_MIN),
    theme = theme(plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
                  plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0.5))) &
  theme(plot.tag = element_text(face = "bold", size = 13))

ggsave(file.path(FIG_DIR, "SuppFig1_QC_panel.pdf"), supp_qc,
       width = 11, height = 8, device = cairo_pdf, dpi = 600)
ggsave(file.path(FIG_DIR, "SuppFig1_QC_panel.png"), supp_qc,
       width = 11, height = 8, dpi = 600)
cat("Supplementary Figure 1 written: SuppFig1_QC_panel.pdf / .png\n")


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
# 5. Compute M-values (asinh-link scale; logit equivalent)
# -----------------------------------------------------------------------------
# Use small offset so beta values at the boundary don't blow up.

eps <- 1e-6
beta_clipped <- pmin(pmax(beta_raw, eps), 1 - eps)
m_raw <- log2(beta_clipped / (1 - beta_clipped))


# -----------------------------------------------------------------------------
# 6. Probe filtering
# -----------------------------------------------------------------------------
# (a) Detection p > 0.01 in any sample
# (b) Cross-reactive probes via Zhou KYCG EPICv2 mask
# (c) Sex chromosomes (will be addressed via covariate later)
# (d) SNP-adjacent (via minfi annotation)
# (e) Invariant: beta > 0.8 or beta < 0.2 in ALL samples

# (a)
keep_detP <- rowSums(detP > 0.01, na.rm = TRUE) == 0
cat("Passing detP filter:", sum(keep_detP), "/", length(keep_detP), "\n")

# (b) Cross-reactive / design-flagged probes from the KYCG EPICv2 mask.
#     This sesameDataGet now resolves via the alternate host (SESAMEDATA_USE_ALT,
#     set in section 1), returning the real curated mask object rather than the
#     cold ExperimentHub copy. The object is a named list of mask sub-DBs (SNP,
#     mapping, non-unique, etc.); unlist+unique collapses them to the full set
#     of design-flagged probe IDs, matching the original conservative filter.
epicv2_mask <- sesameDataGet("KYCG.EPICv2.Mask.20230314")
cross_reactive <- unique(unlist(epicv2_mask, use.names = FALSE))
cat("Cross-reactive probes (KYCG EPICv2 mask):", length(cross_reactive), "\n")

# (c) sex chromosomes
annot_epic <- as.data.frame(
  getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
) %>%
  rownames_to_column("probeID")
sex_chr_probes <- annot_epic$probeID[annot_epic$chr %in% c("chrX", "chrY")]

# (d) SNP-adjacent: use sesame's masked probes (already in epicv2_mask if SNP)
# minfi's dropLociWithSnps is for minfi GenomicRatioSet, which we're not using
# here. The KYCG mask already includes SNP-affected probes.

# (e) Invariant
all_high <- rowSums(beta_raw > 0.8, na.rm = TRUE) == ncol(beta_raw)
all_low  <- rowSums(beta_raw < 0.2, na.rm = TRUE) == ncol(beta_raw)
keep_var <- !(all_high | all_low)
cat("Passing variance/saturation filter:", sum(keep_var), "/", length(keep_var), "\n")

keep_probes <- keep_detP & keep_var &
  !(rownames(beta_raw) %in% cross_reactive) &
  !(rownames(beta_raw) %in% sex_chr_probes)
cat("Probes retained after all filters:", sum(keep_probes), "\n")

beta_filt <- beta_raw[keep_probes, ]
m_filt    <- m_raw[keep_probes, ]


# -----------------------------------------------------------------------------
# 7. Collapse EPIC v2.0 replicate probes to base IDs
# -----------------------------------------------------------------------------
# EPIC v2 includes 5-6k probes designed in duplicate or higher with suffixes
# like _TC11, _BC21. Reference panels use base cg########## IDs, so we
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


# -----------------------------------------------------------------------------
# 8. Cell-type deconvolution (EpiDISH hepidish)
# -----------------------------------------------------------------------------
# Two-stage reference:
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
# working `cell_props` keyed by Sample_Name, since section 9 (raw-PCA
# correlation) and section 10 (removeBatchEffect covariates) match on
# colnames(m_collapsed), which are Sample_Names.
cell_props_out <- cell_props %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID),
            by = "Sample_Name") %>%
  dplyr::select(Subject_ID, everything(), -Sample_Name)

cat("\nDeconvolution results (rounded to 3dp):\n")
print(cell_props_out %>% mutate(across(where(is.numeric), ~round(.x, 3))))
write_csv(cell_props_out, file.path(TBL_DIR, "cell_composition_estimates.csv"))

# -----------------------------------------------------------------------------
# Table 1: cell-type composition collapsed to Epi / Fib / Neutrophil / Other immune
# -----------------------------------------------------------------------------
# Collapse the 9 deconvolved fractions into four readable categories. "Other
# immune" sums every immune lineage except neutrophils.
other_immune_cols <- intersect(c("B", "NK", "CD4T", "CD8T", "Mono", "Eosino"),
                               names(cell_props))
table1_props <- cell_props %>%
  transmute(
    Sample_Name,
    Epithelial     = Epi,
    Fibroblast     = Fib,
    Neutrophil     = Neutro,
    `Other immune` = rowSums(across(all_of(other_immune_cols)))
  ) %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID, Sample_Group),
            by = "Sample_Name") %>%
  arrange(Sample_Group, Subject_ID) %>%
  transmute(`Subject ID`   = Subject_ID,
            Group          = Sample_Group,
            Epithelial     = round(100 * Epithelial, 1),
            Fibroblast     = round(100 * Fibroblast, 1),
            Neutrophil     = round(100 * Neutrophil, 1),
            `Other immune` = round(100 * `Other immune`, 1))

save_pub_table(
  table1_props,
  file.path(TBL_DIR, "Table1_cell_composition_collapsed"),
  caption = "Table 1. Estimated cell-type composition by subject (%).",
  note = paste("Values are EpiDISH (hepidish, RPC) proportion estimates,",
               "expressed as percentages. Other immune = B + NK + CD4+ T +",
               "CD8+ T + monocytes + eosinophils. Rows may not sum to exactly",
               "100 due to rounding."))

# -----------------------------------------------------------------------------
# 9. Diagnostic: how much does composition explain the raw methylome?
# -----------------------------------------------------------------------------
# This step also produces the raw PCA needed to demonstrate that PC1 of the
# unadjusted matrix is dominated by composition before we adjust it out.

probe_var_raw <- rowVars(m_collapsed, na.rm = TRUE)
top5_raw_idx  <- which(probe_var_raw >= quantile(probe_var_raw, 0.95, na.rm = TRUE))
m_raw_top5    <- m_collapsed[top5_raw_idx, ]

pca_raw <- prcomp(t(m_raw_top5), center = TRUE, scale. = FALSE)
var_exp_raw <- (pca_raw$sdev^2) / sum(pca_raw$sdev^2)
cat("Raw-data PC variance explained (first 5):",
    paste0(round(var_exp_raw[1:5] * 100, 1), "%", collapse = ", "), "\n")

cell_cols <- setdiff(names(cell_props), "Sample_Name")
pc1_vec <- pca_raw$x[cell_props$Sample_Name, "PC1"]

cor_table <- map_dfr(cell_cols, function(ct) {
  res <- cor.test(cell_props[[ct]], pc1_vec)
  tibble(cell_type = ct,
         pearson_r = unname(res$estimate),
         p_value   = res$p.value)
}) %>% arrange(desc(abs(pearson_r)))

cat("\nCorrelation of each cell-type proportion with raw-methylome PC1:\n")
print(cor_table %>% mutate(pearson_r = round(pearson_r, 3),
                           p_value   = signif(p_value, 3)))
write_csv(cor_table, file.path(TBL_DIR, "composition_vs_raw_pc1.csv"))

# -----------------------------------------------------------------------------
# Table 2: composition vs unadjusted PC1 (publication-ready)
# -----------------------------------------------------------------------------
ct_full <- c(Epi = "Epithelial", Fib = "Fibroblast", Neutro = "Neutrophil",
             Mono = "Monocyte", Eosino = "Eosinophil", CD8T = "CD8+ T cell",
             CD4T = "CD4+ T cell", NK = "NK cell", B = "B cell")
fmt_p <- function(p) ifelse(p < 0.001,
                            formatC(p, format = "e", digits = 1),
                            formatC(p, format = "f", digits = 3))
table2_pc1 <- cor_table %>%
  transmute(`Cell type` = dplyr::recode(cell_type, !!!ct_full),
            `Pearson r` = round(pearson_r, 3),
            `p-value`   = fmt_p(p_value))

save_pub_table(
  table2_pc1,
  file.path(TBL_DIR, "Table2_composition_vs_rawPC1"),
  caption = paste0("Table 2. Correlation of cell-type proportions with the ",
                   "first principal component of the unadjusted methylome ",
                   "(PC1 = ", round(var_exp_raw[1] * 100, 1), "% of variance)."),
  note = paste("Pearson correlation between EpiDISH-estimated proportions and",
               "PC1 of the top-5% most variable CpGs in the unadjusted",
               "(pre-composition-adjustment) methylome. Cell types are ordered",
               "by absolute correlation."))

# -----------------------------------------------------------------------------
# Figure 1: 2D PCA of the unadjusted (pre-composition-adjustment) methylome
# -----------------------------------------------------------------------------
# PC1 of the raw methylome tracks epithelial content (Table 2), so points are
# colored by group. Subject-ID labels with leader stems are added to match the
# composition-adjusted PCA (section 11).
fig1_df <- as.data.frame(pca_raw$x[, 1:2]) %>%
  rownames_to_column("Sample_Name") %>%
  left_join(dplyr::select(sample_sheet, Sample_Name, Subject_ID, Sample_Group),
            by = "Sample_Name") %>%
  mutate(label = gsub("_", "", Subject_ID))

fig1_pca <- ggplot(fig1_df, aes(PC1, PC2, color = Sample_Group, label = label)) +
  geom_point(size = 3.4) +
  geom_text_repel(size = 2.7, show.legend = FALSE,
                  min.segment.length = 0,        # draw a stem for every label
                  segment.size = 0.3, segment.color = "grey55",
                  box.padding = 0.5, point.padding = 0.25,
                  max.overlaps = Inf, seed = 1) +
  scale_color_manual(values = c(ARDS = "#D55E00", Control = "#56B4E9"), name = "Group") +
  labs(x = paste0("PC1 (", round(var_exp_raw[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_exp_raw[2] * 100, 1), "%)"),
       title = "Unadjusted Methylome PCA") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(file.path(FIG_DIR, "Figure1_PCA_unadjusted.pdf"), fig1_pca,
       width = 7, height = 5, device = cairo_pdf, dpi = 600)
ggsave(file.path(FIG_DIR, "Figure1_PCA_unadjusted.png"), fig1_pca,
       width = 7, height = 5, dpi = 600)
cat("Figure 1 written: Figure1_PCA_unadjusted.pdf / .png\n")


# -----------------------------------------------------------------------------
# 10. Composition-adjusted M-values
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

# (b) drop near-zero-variance columns. Threshold of 1e-4 in variance catches
#     CD8T and Eosino when they are zero or near-zero in all samples.
candidate_cols <- setdiff(cell_cols, drop_col)
col_vars <- apply(prop_mat[, candidate_cols, drop = FALSE], 2, var)
near_zero_var_cols <- names(col_vars)[col_vars < 1e-4]
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
# 11. PCA on adjusted matrix
# -----------------------------------------------------------------------------
# Recompute top 5% variable CpGs on the adjusted matrix (the set will differ
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
  mutate(label = gsub("_", "", Subject_ID))

# Targeted nudge: pull the ARDS026 label up and to the left, clear of its own
# dot in the upper-left corner. min.segment.length = 0 keeps a leader stem from
# the dot to the moved label. Tune the magnitudes (or flip the signs) to taste.
# The panel expansion gives the moved label room so it is not clipped at the edge.
# Per-label manual nudges for THIS layout. ggrepel won't reroute a stem that
# crosses a point, so we move offending labels by hand. Keyed by label; any
# label not listed gets 0. Tune magnitudes / flip signs to taste.
nudge_x_map <- c(ARDS026 = -12, ARDS004 =  16, Control001 = -10, Control022 =  10)
nudge_y_map <- c(ARDS026 =  12, ARDS004 = -10, Control001 =   8, Control022 =  -8)

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
  scale_color_manual(values = c("Control" = "#56B4E9", "ARDS" = "#D55E00")) +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(x = paste0("PC1 (", round(var_exp_adj[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(var_exp_adj[2] * 100, 1), "%)"),
       title = "PCA on composition-adjusted methylome") +
  theme_classic(base_size = 12)

ggsave(file.path(FIG_DIR, "pca_adjusted_2d.pdf"),
       p_pca, width = 7, height = 5, device = cairo_pdf, dpi = 600)

# -----------------------------------------------------------------------------
# Table 3: PCA variance explained, unadjusted vs composition-adjusted
# -----------------------------------------------------------------------------
# Publication-ready version of pca_variance_comparison.csv. Percentages are
# forced to one decimal so the .docx column reads cleanly (e.g. "4.0" rather
# than "4").
fmt1 <- function(x) sprintf("%.1f", x)

pc1_raw_pct <- fmt1(var_exp_raw[1] * 100)
pc1_adj_pct <- fmt1(var_exp_adj[1] * 100)

table3_var <- tibble(
  `Principal component` = paste0("PC", 1:10),
  `Unadjusted (%)`      = fmt1(var_exp_raw[1:10] * 100),
  `Adjusted (%)`        = fmt1(var_exp_adj[1:10] * 100)
)

save_pub_table(
  table3_var,
  file.path(TBL_DIR, "Table3_pca_variance"),
  caption = paste0("Table 3. Variance explained by the first 10 principal ",
                   "components of the unadjusted versus composition-adjusted ",
                   "methylome."),
  note = paste0("Principal components were computed on the top 5% most variable ",
                "CpGs in each matrix. In the unadjusted methylome a single ",
                "component (PC1) accounts for ", pc1_raw_pct, "% of variance and ",
                "tracks epithelial cell fraction (see Table 2); after composition ",
                "adjustment variance is redistributed across many components, ",
                "with PC1 falling to ", pc1_adj_pct, "%."))
# -----------------------------------------------------------------------------
# 12. Hierarchical clustering with silhouette-based k selection
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

# Dendrogram
cairo_pdf(file.path(FIG_DIR, "dendrogram_adjusted.pdf"), width = 10, height = 5)
plot(sample_hc_adj,
     labels = sample_sheet$Subject_ID[match(sample_hc_adj$labels,
                                            sample_sheet$Sample_Name)],
     main = paste0("Composition-adjusted clustering (k = ", chosen_k, ")"),
     xlab = "", sub = "", cex = 0.8)
dev.off()



# -----------------------------------------------------------------------------
# Supplementary Table: silhouette width by k (cluster-structure diagnostic)
# -----------------------------------------------------------------------------
# Publication-ready version of silhouette_by_k.csv.
suppTable_sil <- sil_table %>%
  transmute(
    `Number of clusters (k)` = k,
    `Mean silhouette width`  = sprintf("%.3f", mean_silhouette))

save_pub_table(
  suppTable_sil,
  file.path(TBL_DIR, "SuppTable_silhouette_by_k"),
  caption = paste0("Supplementary Table. Mean silhouette width for k = 2 to 6 ",
                   "clusters on the composition-adjusted methylome."),
  note = paste0("Silhouette widths were computed on Ward.D2 hierarchical ",
                "clustering of the top 5% most variable composition-adjusted ",
                "CpGs. Interpretation thresholds: > 0.5 strong, 0.25 to 0.5 ",
                "weak, < 0.25 essentially absent. The best-supported solution ",
                "(k = ", best_k, ", width ", sprintf("%.3f", best_sil), ") ",
                "falls below 0.25, indicating that no robust discrete subgroup ",
                "structure survives adjustment for cell composition."))


# -----------------------------------------------------------------------------
# Annotated dendrogram: cluster assignments (bottom strip)
# -----------------------------------------------------------------------------
# Honest depiction of the (weak) k = chosen_k partition. The tree is left BLACK
# so branch color does not overstate the separation; the k = chosen_k cut is
# drawn; and a single strip beneath the leaves encodes the cluster assignment
# stored in Adjusted_Group. Diagnosis is readable from the Subject_ID labels.
# Pairs with SuppTable_silhouette_by_k (width 0.233 at k = 2, below the 0.25
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

# (c) Palette (colorblind-safe; deliberately not the diagnosis colors)
clust_palette <- c("#009E73", "#CC79A7", "#0072B2", "#E69F00", "#000000", "#999999")
clust_cols    <- clust_palette[seq_len(chosen_k)]

# (d) Tree stays black; only set label size. The strip carries the assignment.
dend <- set(dend, "labels_cex", 0.7)

# (e) Relabel leaves to Subject_ID for display (leaf order unchanged)
labels(dend) <- unname(setNames(sample_sheet$Subject_ID,
                                sample_sheet$Sample_Name)[labels(dend)])

# (f) Strip colors in leaf order, from the stored assignment
strip_cols <- clust_cols[clusters_leaf]

# (g) Cut height for the k-cluster solution (midway between bracketing merges)
n_leaves <- length(leaf_samples)
hts      <- sort(sample_hc_adj$height)
cut_h    <- mean(c(hts[n_leaves - chosen_k], hts[n_leaves - chosen_k + 1]))

# (h) Draw helper, reused for PDF and PNG
draw_dendro <- function() {
  par(mar = c(9, 4.5, 5, 2), xpd = NA)   # taller top margin; xpd lets legend sit above
  plot(dend,
       main = paste0("Composition-adjusted clustering (k = ", chosen_k, ")"),
       ylab = "Height", xlab = "", sub = "")
  abline(h = cut_h, lty = 2, col = "grey40")
  colored_bars(colors = strip_cols, dend = dend, rowLabels = "Cluster",
               sort_by_labels_order = FALSE)
  legend("topright", inset = c(0, -0.12), bty = "n", cex = 0.8,
         title = paste0("Cluster (k = ", chosen_k, ")"),
         legend = paste0("Group ", seq_len(chosen_k)),
         fill   = clust_cols)
}

cairo_pdf(file.path(FIG_DIR, "dendrogram_adjusted_annotated.pdf"),
          width = 11, height = 6)
draw_dendro(); dev.off()
png(file.path(FIG_DIR, "dendrogram_adjusted_annotated.png"),
    width = 11, height = 6, units = "in", res = 600)
draw_dendro(); dev.off()
cat("Annotated dendrogram written: dendrogram_adjusted_annotated.pdf / .png\n")

# Checkpoint: everything needed to resume from section 13 onward without
# re-running openSesame / hepidish / removeBatchEffect.
saveRDS(list(
  sample_sheet = sample_sheet, cell_props = cell_props, cell_cols = cell_cols,
  prop_for_adj = prop_for_adj, m_collapsed = m_collapsed,
  beta_collapsed = beta_collapsed, m_adjusted = m_adjusted,
  m_adj_top5 = m_adj_top5, pca_adj = pca_adj, var_exp_adj = var_exp_adj,
  var_exp_raw = var_exp_raw, sample_hc_adj = sample_hc_adj,
  chosen_k = chosen_k, annot_epic = annot_epic
), file.path(RDS_DIR, "checkpoint_post_clustering.rds"))

# -----------------------------------------------------------------------------
# 13. Clinical outcome analysis by adjusted groups
# -----------------------------------------------------------------------------
# Compare clinical variables across adjusted clusters. Continuous variables
# by Wilcoxon (Kruskal-Wallis if > 2 groups); categorical by Fisher exact.
#
# Available clinical columns (from external metadata, joined in section 2):
#   Subject, Group, Age, Sex, Race, PARDS_severity, PELOD, VFD,
#   Principal_Comorbidity, Acute_Condition, COVID, Direct_Lung_Injury,
#   Infectious_Agent, Outcome, Methyl_Group
#
# Some derived variables are useful for testing:
#   - immunocompromised flag derived from Principal_Comorbidity
#   - viral_present flag derived from Infectious_Agent
#   - any_PARDS flag (collapse severity to binary) for sanity check

# Derive convenience flags if the source columns exist
if ("Principal_Comorbidity" %in% names(sample_sheet)) {
  sample_sheet <- sample_sheet %>%
    mutate(Immunocompromised_derived = grepl("immun|onc|transplant|hsct|chemo",
                                             Principal_Comorbidity,
                                             ignore.case = TRUE))
}
if ("Infectious_Agent" %in% names(sample_sheet)) {
  sample_sheet <- sample_sheet %>%
    mutate(Viral_present_derived = grepl("virus|viral|rsv|rhino|adeno|covid|influenza|metapneumo|parainfluenza",
                                         Infectious_Agent,
                                         ignore.case = TRUE))
}

continuous_vars  <- c("Age", "VFD", "PELOD")
categorical_vars <- c("Sex", "Race", "PARDS_severity", "Outcome",
                      "Direct_Lung_Injury", "COVID",
                      "Principal_Comorbidity", "Acute_Condition",
                      "Infectious_Agent",
                      "Immunocompromised_derived", "Viral_present_derived",
                      "Group", "Methyl_Group")

cont_avail <- intersect(continuous_vars,  names(sample_sheet))
cat_avail  <- intersect(categorical_vars, names(sample_sheet))
cat("\nContinuous variables for testing:",
    if (length(cont_avail)) paste(cont_avail, collapse = ", ") else "(none)", "\n")
cat("Categorical variables for testing:",
    if (length(cat_avail))  paste(cat_avail,  collapse = ", ") else "(none)", "\n")

compare_continuous <- function(df, var, grp) {
  v <- df[[var]]
  g <- df[[grp]]
  if (sum(!is.na(v)) < 4 || length(unique(na.omit(g))) < 2) {
    return(tibble(variable = var, test = "skipped", p_value = NA_real_))
  }
  
  test_result <- tryCatch({
    if (length(unique(na.omit(g))) == 2) {
      list(t = wilcox.test(v ~ g, exact = FALSE), name = "Mann-Whitney U")
    } else {
      list(t = kruskal.test(v ~ as.factor(g)), name = "Kruskal-Wallis")
    }
  }, error = function(e) {
    list(t = NULL, name = paste("failed:", conditionMessage(e)))
  }, warning = function(w) {
    # Some tests issue ties warnings; still want the result
    list(t = suppressWarnings(
      if (length(unique(na.omit(g))) == 2) wilcox.test(v ~ g, exact = FALSE)
      else kruskal.test(v ~ as.factor(g))),
      name = if (length(unique(na.omit(g))) == 2) "Mann-Whitney U" else "Kruskal-Wallis")
  })
  
  tibble(
    variable = var,
    test     = test_result$name,
    p_value  = if (!is.null(test_result$t)) test_result$t$p.value else NA_real_
  )
}

compare_categorical <- function(df, var, grp) {
  tbl <- table(df[[var]], df[[grp]])
  if (any(dim(tbl) < 2)) {
    return(tibble(variable = var, test = "skipped", p_value = NA_real_))
  }
  
  test_result <- tryCatch({
    fisher.test(tbl, simulate.p.value = TRUE, B = 1e5)
  }, error = function(e) {
    NULL
  })
  
  tibble(
    variable = var,
    test     = if (!is.null(test_result)) "Fisher exact" else "failed",
    p_value  = if (!is.null(test_result)) test_result$p.value else NA_real_
  )
}
cont_results <- if (length(cont_avail)) {
  map_dfr(cont_avail, ~compare_continuous(sample_sheet, .x, "Adjusted_Group"))
} else tibble()
cat_results <- if (length(cat_avail)) {
  map_dfr(cat_avail, ~compare_categorical(sample_sheet, .x, "Adjusted_Group"))
} else tibble()

clinical_results <- bind_rows(cont_results, cat_results) %>%
  arrange(p_value)
cat("\nClinical comparisons by Adjusted_Group:\n")
print(clinical_results %>% mutate(p_value = signif(p_value, 3)))
write_csv(clinical_results,
          file.path(TBL_DIR, "clinical_comparisons_by_adjusted_group.csv"))

# -----------------------------------------------------------------------------
# Table 4: characteristics by adjusted group (publication-ready)
# -----------------------------------------------------------------------------
# Continuous: median [Q1, Q3]; categorical: n (column %). Last column is the
# p-value from the SAME test section 13 used (Mann-Whitney / Kruskal-Wallis for
# continuous; simulated Fisher exact for categorical). One p-value per variable,
# placed on the variable's header row for categoricals.
set.seed(20260101)   # reproducible simulated Fisher p-values

grp_var    <- "Adjusted_Group"
grp_levels <- sort(unique(na.omit(sample_sheet[[grp_var]])))
grp_n      <- vapply(grp_levels,
                     function(g) sum(sample_sheet[[grp_var]] == g, na.rm = TRUE),
                     integer(1))
# Match the figure legend ("Group 1/2") in the column headers
grp_disp    <- sub("^AdjGroup", "Group ", grp_levels)
grp_headers <- setNames(paste0(grp_disp, " (n=", grp_n, ")"), grp_levels)

ordered_levels <- list(
  `PARDS severity` = c("None", "Mild", "Moderate", "Severe")
)

fmt_pval <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) return("<0.001")
  formatC(p, format = "f", digits = 3)
}
fmt_med_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return("-")
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE, type = 7)
  sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3])
}
cont_p <- function(v, g) {
  ok <- !is.na(v) & !is.na(g); v <- v[ok]; g <- droplevels(as.factor(g[ok]))
  if (length(v) < 4 || nlevels(g) < 2) return(NA_real_)
  suppressWarnings(if (nlevels(g) == 2) wilcox.test(v ~ g, exact = FALSE)$p.value
                   else kruskal.test(v ~ g)$p.value)
}
cat_p <- function(x, g) {
  tbl <- table(x, g)
  if (any(dim(tbl) < 2)) return(NA_real_)
  suppressWarnings(fisher.test(tbl, simulate.p.value = TRUE, B = 1e5)$p.value)
}
# Constant column order for EVERY row, so bind_rows has nothing to misalign.
row_cols <- c("Characteristic", unname(grp_headers[grp_levels]), "p-value")

make_row <- function(characteristic, group_vals, pval_str = "") {
  vals <- c(characteristic,
            as.character(group_vals[grp_levels]),   # group cols, in grp order
            pval_str)
  as_tibble(setNames(as.list(vals), row_cols))      # name the list first, THEN tibble
}
cont_row <- function(var, label) {
  v <- suppressWarnings(as.numeric(sample_sheet[[var]])); g <- sample_sheet[[grp_var]]
  gv <- setNames(vapply(grp_levels, function(L) fmt_med_iqr(v[g == L & !is.na(g)]),
                        character(1)), grp_levels)
  make_row(paste0(label, ", median [IQR]"), gv, fmt_pval(cont_p(v, g)))
}
cat_rows <- function(var, label) {
  g <- sample_sheet[[grp_var]]; x <- as.character(sample_sheet[[var]])
  x <- dplyr::recode(x,
                     "TRUE" = "Yes", "FALSE" = "No",
                     "Y" = "Yes", "N" = "No",
                     "yes" = "Yes", "no" = "No",
                     .default = x)
  keep <- !is.na(x) & !is.na(g); x <- x[keep]; gg <- g[keep]
  if (!length(x)) return(tibble())
  p   <- if (length(unique(x)) >= 2) cat_p(x, gg) else NA_real_
  blank <- setNames(rep("", length(grp_levels)), grp_levels)
  hdr <- make_row(label, blank, fmt_pval(p))
  
  # Order levels: use the explicit order for this label if defined, keeping only
  # levels actually present, then append any unlisted levels alphabetically.
  present <- unique(x)
  ord     <- ordered_levels[[label]]
  level_order <- if (is.null(ord)) sort(present)
  else c(intersect(ord, present), sort(setdiff(present, ord)))
  
  levs <- map_dfr(level_order, function(Lv) {
    gv <- setNames(vapply(grp_levels, function(L) {
      d <- sum(gg == L); n <- sum(gg == L & x == Lv)
      if (d == 0) "-" else sprintf("%d (%.0f%%)", n, 100 * n / d)
    }, character(1)), grp_levels)
    make_row(paste0("\u2003", Lv), gv, "")
  })
  bind_rows(hdr, levs)
}

labels_map <- c(Age = "Age (years)", VFD = "Ventilator-free days",
                PELOD = "PELOD score", Sex = "Sex", Race = "Race",
                PARDS_severity = "PARDS severity", Outcome = "Outcome",
                Direct_Lung_Injury = "Direct lung injury", COVID = "COVID-19",
                Immunocompromised_derived = "Immunocompromised",
                Viral_present_derived = "Viral pathogen detected",
                Sample_Group = "Diagnosis (ARDS / Control)",
                Methyl_Group = "Original methylation group")
lab <- function(v) ifelse(is.na(labels_map[v]), v, labels_map[v])

cont_vars <- intersect(c("Age", "VFD", "PELOD"), names(sample_sheet))
cat_vars  <- intersect(c("Sex", "Race", "PARDS_severity", "Outcome",
                         "Direct_Lung_Injury", "COVID",
                         "Immunocompromised_derived", "Viral_present_derived",
                         "Sample_Group"), names(sample_sheet))

table4_char <- bind_rows(
  map_dfr(cont_vars, ~cont_row(.x, lab(.x))),
  map_dfr(cat_vars,  ~cat_rows(.x, lab(.x))))

save_pub_table(
  table4_char,
  file.path(TBL_DIR, "Table4_characteristics_by_adjusted_group"),
  caption = "Table 4. Clinical characteristics by composition-adjusted cluster.",
  note = paste0("Continuous variables are median [Q1, Q3]; categorical ",
                "variables are n (column %). p-values compare adjusted groups ",
                "by Mann-Whitney U or Kruskal-Wallis (continuous) and simulated ",
                "Fisher exact test (categorical), and are descriptive and ",
                "unadjusted for multiple comparisons. Groups are derived from ",
                "hierarchical clustering of the composition-adjusted methylome."))

# -----------------------------------------------------------------------------
# 13b. Cell-type proportions vs clinical variables and outcomes
# -----------------------------------------------------------------------------
# Place this AFTER section 13 (it reuses the Immunocompromised_derived and
# Viral_present_derived flags created there) and BEFORE section 14.
#
# This asks a question that is distinct from the methylation analysis: now that
# composition has been estimated, do the *cell-type proportions themselves*
# differ by phenotype (ARDS vs control, immunocompromised vs not, ...) or track
# with continuous outcomes (VFD, PELOD)?
#
# Caveats baked into the design:
#   - Proportions are compositional (sum to 1), so these are marginal,
#     per-cell-type tests, not a joint compositional model. With small n and
#     several structurally-zero cell types (CD8T, Eosino, often Fib), a CLR or
#     Dirichlet model is unstable, so we test each cell type marginally and
#     control the FDR across tests.
#   - Cell types detected in too few samples carry no information and are
#     dropped before testing (same spirit as the near-zero-variance filter in
#     section 10). Otherwise a cell type that is 0 in 23/24 samples generates a
#     meaningless "significant" hit off a single nonzero value.
#   - Continuous outcomes: Spearman. Two-group: Mann-Whitney. >2 groups:
#     Kruskal-Wallis.

# Merge proportions with clinical data (cell_props is keyed by Sample_Name)
prop_clin <- cell_props %>%
  left_join(sample_sheet, by = "Sample_Name")

cell_cols <- setdiff(names(cell_props), "Sample_Name")

# Drop low-prevalence / near-constant cell types before testing
prop_only    <- as.matrix(cell_props[, cell_cols, drop = FALSE])
frac_nonzero <- colMeans(prop_only > 0, na.rm = TRUE)
ct_var       <- apply(prop_only, 2, var, na.rm = TRUE)
testable_ct  <- cell_cols[frac_nonzero >= 0.20 & ct_var >= 1e-4]
dropped_ct   <- setdiff(cell_cols, testable_ct)
if (length(dropped_ct) > 0) {
  cat("Cell types excluded (low prevalence / near-constant):",
      paste(dropped_ct, collapse = ", "), "\n")
}
cat("Cell types tested:", paste(testable_ct, collapse = ", "), "\n")

# Outcome / grouping variables (reuse section 13 definitions where present)
cont_outcomes <- intersect(c("Age", "VFD", "PELOD", "Total_Days_Intubated"),
                           names(prop_clin))
cat_outcomes  <- intersect(c("Sample_Group", "Group", "Outcome",
                             "PARDS_severity", "Sex", "Direct_Lung_Injury",
                             "COVID", "Immunocompromised_derived",
                             "Viral_present_derived", "Methyl_Group"),
                           names(prop_clin))
cat("Continuous outcomes:",
    if (length(cont_outcomes)) paste(cont_outcomes, collapse = ", ") else "(none)", "\n")
cat("Grouping variables:  ",
    if (length(cat_outcomes))  paste(cat_outcomes,  collapse = ", ") else "(none)", "\n")

# --- continuous outcomes: Spearman correlation ---
cont_assoc <- map_dfr(testable_ct, function(ct) {
  map_dfr(cont_outcomes, function(out) {
    x <- prop_clin[[ct]]; y <- suppressWarnings(as.numeric(prop_clin[[out]]))
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 4) return(tibble())
    r <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman",
                                   exact = FALSE))
    tibble(cell_type = ct, variable = out, test = "Spearman",
           n = sum(ok), estimate = unname(r$estimate),
           p_value = r$p.value, note = "rho")
  })
})

# --- categorical / grouping variables: Mann-Whitney or Kruskal-Wallis ---
cat_assoc <- map_dfr(testable_ct, function(ct) {
  map_dfr(cat_outcomes, function(grp) {
    x <- prop_clin[[ct]]; g <- prop_clin[[grp]]
    ok <- is.finite(x) & !is.na(g)
    g  <- droplevels(as.factor(g[ok])); x <- x[ok]
    if (nlevels(g) < 2 || length(x) < 4) return(tibble())
    if (nlevels(g) == 2) {
      w    <- suppressWarnings(wilcox.test(x ~ g, exact = FALSE))
      meds <- tapply(x, g, median, na.rm = TRUE)
      tibble(cell_type = ct, variable = grp, test = "Mann-Whitney",
             n = length(x), estimate = unname(meds[2] - meds[1]),
             p_value = w$p.value,
             note = paste0("median higher in ", names(meds)[which.max(meds)]))
    } else {
      k <- kruskal.test(x ~ g)
      tibble(cell_type = ct, variable = grp, test = "Kruskal-Wallis",
             n = length(x), estimate = NA_real_,
             p_value = k$p.value, note = paste0(nlevels(g), " groups"))
    }
  })
})

# FDR control across all proportion tests. To correct within each clinical
# variable instead, replace the mutate with a group_by(variable) first.
prop_assoc <- bind_rows(cont_assoc, cat_assoc) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

cat("\nCell-type proportion vs clinical/outcome associations (BH across all tests):\n")
print(prop_assoc %>%
        mutate(estimate = round(estimate, 3),
               p_value  = signif(p_value, 3),
               p_adj_BH = signif(p_adj_BH, 3)),
      n = Inf)
write_csv(prop_assoc, file.path(TBL_DIR, "cellprop_outcome_associations.csv"))

# --- visualize proportions across the headline grouping variables ---
plot_groups <- intersect(c("Sample_Group", "Immunocompromised_derived",
                           "Outcome"), names(prop_clin))

prop_long <- prop_clin %>%
  dplyr::select(Sample_Name, all_of(testable_ct), all_of(plot_groups)) %>%
  pivot_longer(all_of(testable_ct), names_to = "cell_type",
               values_to = "proportion")

for (grp in plot_groups) {
  df <- prop_long %>% filter(!is.na(.data[[grp]]))
  if (length(unique(df[[grp]])) < 2) next
  p <- ggplot(df, aes(x = factor(.data[[grp]]), y = proportion,
                      fill = factor(.data[[grp]]))) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
    facet_wrap(~ cell_type, scales = "free_y") +
    labs(x = grp, y = "Estimated proportion",
         title = paste0("Cell-type proportions by ", grp)) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(FIG_DIR, paste0("cellprop_by_", grp, ".pdf")),
         p, width = 8, height = 6, device = cairo_pdf)
}


# -----------------------------------------------------------------------------
# 14. Genomic annotation and enrichment setup
# -----------------------------------------------------------------------------
# Builds the promoter/body/intergenic probe annotation, the MSigDB gene-set
# collection, and the run_ora() helper used by the section 15 enrichment.

# Reuse probe annotation from earlier
probe_anno_v2 <- annot_epic %>%
  mutate(base_id = sub("_[A-Z0-9]+$", "", probeID)) %>%
  group_by(base_id) %>% slice_head(n = 1) %>% ungroup() %>%
  rename(probeID_full = probeID, probeID = base_id) %>%
  filter(probeID %in% rownames(m_collapsed))

# Build a region annotation: promoter (-1500 / +200), body, intergenic
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- suppressMessages(genes(txdb))
genes_gr$symbol <- mapIds(org.Hs.eg.db, keys = genes_gr$gene_id,
                          column = "SYMBOL", keytype = "ENTREZID",
                          multiVals = "first")
standard_chrs <- paste0("chr", c(1:22, "X", "Y"))
genes_gr <- genes_gr[seqnames(genes_gr) %in% standard_chrs]
seqlevels(genes_gr, pruning.mode = "coarse") <- standard_chrs

promoter_gr <- trim(promoters(genes_gr, upstream = 1500, downstream = 200))

body_gr <- genes_gr[width(genes_gr) > 200]
is_plus  <- as.character(strand(body_gr)) == "+"
ns <- start(body_gr); ne <- end(body_gr)
ns[is_plus]  <- ns[is_plus] + 200
ne[!is_plus & as.character(strand(body_gr)) == "-"] <-
  ne[!is_plus & as.character(strand(body_gr)) == "-"] - 200
ranges(body_gr) <- IRanges(start = ns, end = ne)
body_gr <- trim(body_gr)

probe_gr <- GRanges(seqnames = probe_anno_v2$chr,
                    ranges = IRanges(probe_anno_v2$pos, probe_anno_v2$pos),
                    probeID = probe_anno_v2$probeID)

prom_hits <- findOverlaps(probe_gr, promoter_gr)
body_hits <- findOverlaps(probe_gr, body_gr)

probe_anno_final <- tibble(probeID = probe_anno_v2$probeID,
                           region = "intergenic", gene = NA_character_)
prom_df <- tibble(probeID = probe_gr$probeID[queryHits(prom_hits)],
                  gene = promoter_gr$symbol[subjectHits(prom_hits)]) %>%
  distinct(probeID, .keep_all = TRUE) %>% mutate(region = "promoter")
body_df <- tibble(probeID = probe_gr$probeID[queryHits(body_hits)],
                  gene = body_gr$symbol[subjectHits(body_hits)]) %>%
  distinct(probeID, .keep_all = TRUE) %>%
  filter(!(probeID %in% prom_df$probeID)) %>% mutate(region = "body")
probe_anno_final <- bind_rows(prom_df, body_df) %>%
  bind_rows(probe_anno_final %>%
              filter(!probeID %in% c(prom_df$probeID, body_df$probeID))) %>%
  distinct(probeID, .keep_all = TRUE)

# Load MSigDB collections (Hallmark, C2 CP, C5 GO)
msig_all <- bind_rows(
  msigdbr(species = "Homo sapiens", category = "H"),
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP"),
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"),
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY"),
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"),
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF"),
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC")
) %>% dplyr::select(gs_name, gene_symbol) %>% distinct()

universe_genes <- unique(na.omit(probe_anno_final$gene))

run_ora <- function(gene_list, label) {
  if (length(gene_list) < 5) return(NULL)
  enricher(gene = gene_list, universe = universe_genes,
           TERM2GENE = msig_all, pAdjustMethod = "BH",
           pvalueCutoff = 0.05, qvalueCutoff = 0.05,
           minGSSize = 10, maxGSSize = 500)
}

# -----------------------------------------------------------------------------
# 15. Biology of adjusted groups: limma DMC analysis
# -----------------------------------------------------------------------------
# DMCs between adjusted groups, fit on collapsed M-values with cell
# proportions as covariates in the same model. This is the recommended
# approach over residualize-then-test.

if (length(unique(sample_sheet$Adjusted_Group)) >= 2) {
  
  grp_factor <- factor(sample_sheet$Adjusted_Group)
  design_dm  <- model.matrix(~ grp_factor + prop_for_adj)
  
  fit  <- lmFit(m_collapsed[, sample_sheet$Sample_Name], design_dm)
  fit2 <- eBayes(fit)
  
  # Test the first non-intercept group factor coefficient (Group2 vs Group1)
  coef_name <- colnames(design_dm)[2]
  cat("\nTesting coefficient:", coef_name, "\n")
  
  dmc_res <- topTable(fit2, coef = coef_name, number = Inf,
                      adjust.method = "BH", sort.by = "none") %>%
    rownames_to_column("probeID")
  
  # Compute delta beta (between the first two groups for interpretability)
  grp_levels <- levels(grp_factor)
  s1 <- sample_sheet$Sample_Name[sample_sheet$Adjusted_Group == grp_levels[1]]
  s2 <- sample_sheet$Sample_Name[sample_sheet$Adjusted_Group == grp_levels[2]]
  delta_beta <- rowMeans(beta_collapsed[, s2, drop = FALSE], na.rm = TRUE) -
    rowMeans(beta_collapsed[, s1, drop = FALSE], na.rm = TRUE)
  dmc_res <- dmc_res %>%
    mutate(delta_beta  = delta_beta[probeID],
           significant = adj.P.Val < 0.05 & abs(delta_beta) >= 0.1,
           direction   = case_when(
             significant & delta_beta >=  0.1 ~ paste0("hyper_in_", grp_levels[2]),
             significant & delta_beta <= -0.1 ~ paste0("hypo_in_",  grp_levels[2]),
             TRUE                             ~ "ns"
           )) %>%
    left_join(probe_anno_final, by = "probeID")
  
  cat("\nDMC summary (adjusted-group contrast, adjusted for cell composition):\n")
  print(table(dmc_res$direction))
  write_csv(dmc_res, file.path(TBL_DIR, "dmc_adjusted_groups.csv"))
  
  # --- Publication count table: how many DMCs between adjusted groups ---
  n_hyper <- sum(dmc_res$direction == paste0("hyper_in_", grp_levels[2]), na.rm = TRUE)
  n_hypo  <- sum(dmc_res$direction == paste0("hypo_in_",  grp_levels[2]), na.rm = TRUE)
  n_sig   <- n_hyper + n_hypo
  cat(sprintf("Significant DMCs (FDR<0.05 & |delta beta|>=0.1): %d (hyper %d / hypo %d)\n",
              n_sig, n_hyper, n_hypo))
  
  grp_disp_dm <- sub("^AdjGroup", "Group ", grp_levels)
  dmc_count_tbl <- tibble(
    Direction  = c(paste0("Hypermethylated in ", grp_disp_dm[2]),
                   paste0("Hypomethylated in ",  grp_disp_dm[2]),
                   "Total significant"),
    `CpGs (n)` = c(n_hyper, n_hypo, n_sig))
  save_pub_table(
    dmc_count_tbl,
    file.path(TBL_DIR, "Table5_dmc_counts_adjusted_groups"),
    caption = paste0("Table 5. Differentially methylated CpGs between ",
                     "composition-adjusted ", grp_disp_dm[1], " and ",
                     grp_disp_dm[2], "."),
    note = paste0("CpGs with Benjamini-Hochberg adjusted p < 0.05 and ",
                  "|delta beta| >= 0.1 from a limma model of M-values on ",
                  "adjusted group plus cell-type composition covariates. ",
                  "Direction is relative to ", grp_disp_dm[2], ". With ",
                  length(s2), " subjects in ", grp_disp_dm[2], ", these counts ",
                  "are exploratory and should be interpreted with caution."))
  
  # ORA for promoter hyper/hypo (diagnostic-verbose; collects results for figure)
  ora_dir_list <- list()
  for (dir in setdiff(unique(dmc_res$direction), "ns")) {
    genes_in_dir <- dmc_res %>%
      filter(significant, region == "promoter", direction == dir, !is.na(gene)) %>%
      pull(gene) %>% unique()
    cat(sprintf("\n[ORA] %s: %d promoter-DMC genes into enrichment\n",
                dir, length(genes_in_dir)))
    ora <- run_ora(genes_in_dir, dir)
    if (is.null(ora)) {
      cat("  -> fewer than 5 input genes; enrichment not run.\n"); next
    }
    ora_df <- as.data.frame(ora)
    if (nrow(ora_df) == 0) {
      cat("  -> enrichment ran but no terms passed p.adjust / q < 0.05.\n"); next
    }
    ora_dir_list[[dir]] <- ora_df
    write_csv(ora_df, file.path(TBL_DIR, paste0("ora_dmc_promoter_", dir, ".csv")))
    cat("Top 10 promoter-", dir, "pathways:\n")
    print(ora_df %>% arrange(p.adjust) %>%
            dplyr::select(Description, GeneRatio, p.adjust, Count) %>% head(10))
  }

  # --- Figure: promoter DMC pathway enrichment dotplot ---
  if (length(ora_dir_list) > 0) {
    tidy_term <- function(x) {
      x <- sub("^(GOBP|GOCC|GOMF|GO|KEGG|REACTOME|HALLMARK|WP|PID|BIOCARTA|NABA)_",
               "", x)
      x <- tolower(gsub("_", " ", x))
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      ifelse(nchar(x) > 45, paste0(substr(x, 1, 42), "..."), x)
    }
    parse_ratio <- function(gr) vapply(strsplit(gr, "/"),
                                       function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
    dir_pretty <- function(d)
      sub("^hypo_in_AdjGroup", "Hypomethylated in Group ",
          sub("^hyper_in_AdjGroup", "Hypermethylated in Group ", d))
    
    TOP_N <- 10
    make_panel <- function(df, dir) {
      d <- df %>% arrange(p.adjust) %>% head(TOP_N) %>%
        transmute(Term = tidy_term(Description),
                  GeneRatio = parse_ratio(GeneRatio),
                  Count, p.adjust) %>%
        mutate(Term = factor(Term, levels = Term[order(GeneRatio)]))
      ggplot(d, aes(GeneRatio, Term, size = Count, color = p.adjust)) +
        geom_point() +
        scale_color_gradient(low = "blue", high = "red",
                             name = "Adjusted p-value") +
        scale_size_continuous(name = "Genes", range = c(2, 7)) +
        scale_y_discrete(expand = expansion(add = 0.6)) +   # room so top/bottom dots aren't clipped
        labs(x = "Gene ratio", y = NULL, title = dir_pretty(dir)) +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 11))
    }
    
    panels <- purrr::imap(ora_dir_list, make_panel)
    p_ora <- patchwork::wrap_plots(panels, ncol = length(panels)) +
      patchwork::plot_annotation(
        title = "Promoter DMC pathway enrichment (composition-adjusted groups)",
        theme = theme(plot.title = element_text(face = "bold", size = 12)))
    
    n_panels <- length(panels)
    ggsave(file.path(FIG_DIR, "Figure_ora_dmc_promoter.pdf"), p_ora,
           width = 7 * n_panels, height = 5.5, device = cairo_pdf, dpi = 600)
    ggsave(file.path(FIG_DIR, "Figure_ora_dmc_promoter.png"), p_ora,
           width = 7 * n_panels, height = 5.5, dpi = 600)
    cat("Pathway enrichment figure written: Figure_ora_dmc_promoter.pdf / .png\n")
  } else {
    cat("No enrichment results to plot.\n")
  }}

# -----------------------------------------------------------------------------
# 15b. Cell-type-specific differential methylation (CellDMC)
# -----------------------------------------------------------------------------
# Place AFTER section 15 (it reuses probe_anno_final from section 14 and the
# Immunocompromised_derived flag from section 13). Uses EpiDISH::CellDMC, which
# is already loaded.
#
# Question answered: within a given cell type (we focus on Epi and Neutro), is
# the methylome different between phenotype groups, holding composition fixed?
# CellDMC fits, at each CpG, an interaction model
#     beta ~ sum_c [ frac_c * (mu_c + delta_c * phenotype) ]
# so delta_c is the cell-type-specific effect of the phenotype. Input is the
# BETA matrix (0-1), NOT M-values.
#
# Power caveats (read before interpreting):
#   - Interaction tests are far hungrier for samples and for *variation in the
#     fractions* than marginal tests. At n~24 only the abundant, variable
#     compartments (Epi, Neutro) realistically support inference; rare types do
#     not. Treat everything here as exploratory / hypothesis-generating.
#   - Estimate sign convention: phenotype is coded 0/1 below, so a positive
#     Estimate means higher methylation in the "1" group within that cell type.

## ---- Cell-fraction matrix for CellDMC (3-compartment: Epi / Neutro / Rest) --
## Seven cell types make the interaction design rank-deficient at this n,
## especially in the smaller phenotype group.
## Collapse all non-target cells into one "Rest" compartment: the mixture still
## sums to 1, so the other cells are accounted for, but the model is full-rank.
frac_full <- as.matrix(
  cell_props[match(colnames(beta_collapsed), cell_props$Sample_Name),
             setdiff(names(cell_props), "Sample_Name"), drop = FALSE]
)
rownames(frac_full) <- colnames(beta_collapsed)

Epi    <- frac_full[, "Epi"]
Neutro <- frac_full[, "Neutro"]
Rest   <- pmax(0, 1 - Epi - Neutro)          # everything else, pooled
frac_cdmc <- cbind(Epi = Epi, Neutro = Neutro, Rest = Rest)
frac_cdmc <- frac_cdmc / rowSums(frac_cdmc)  # guard tiny rounding drift

report_ct <- c("Epi", "Neutro")              # Rest is a nuisance compartment
## ---- Define contrasts -------------------------------------------------------
## Each function returns a 0/1 vector aligned to rownames(frac_cdmc); NA marks
## samples excluded from that contrast. EDIT the level/regex strings to match
## your actual sample_sheet coding.
ss <- sample_sheet[match(rownames(frac_cdmc), sample_sheet$Sample_Name), ]

contrasts_list <- list(
  ARDS_vs_Control = function(s) {
    ifelse(s$Sample_Group == "ARDS", 1L,
           ifelse(s$Sample_Group == "Control", 0L, NA_integer_))
  },
  Nonsurvivor_vs_Survivor = function(s) {       # 1 = non-survivor
    out <- tolower(as.character(s$Outcome))
    ifelse(grepl("died|death|non.?surviv|expired|mortal", out), 1L,
           ifelse(grepl("surviv|alive|discharg", out),       0L, NA_integer_))
  },
  IC_vs_nonIC = function(s) {                    # 1 = immunocompromised
    if (!"Immunocompromised_derived" %in% names(s))
      return(rep(NA_integer_, nrow(s)))
    as.integer(s$Immunocompromised_derived)
  }
)

## ---- Run CellDMC per contrast ----------------------------------------------
celldmc_summary <- list()

for (cname in names(contrasts_list)) {
  pheno <- contrasts_list[[cname]](ss)
  ok    <- !is.na(pheno)
  tab   <- table(pheno[ok])
  cat("\n=== CellDMC contrast:", cname, "===\n")
  cat("Group sizes (0/1):",
      paste(names(tab), tab, sep = "=", collapse = ", "), "\n")
  
  if (length(tab) < 2 || min(tab) < 4) {
    cat("Skipping: fewer than 4 samples in a group.\n"); next
  }
  
  beta_sub <- beta_collapsed[, ok, drop = FALSE]
  frac_sub <- frac_cdmc[ok, , drop = FALSE]
  pheno_v  <- pheno[ok]
  
  # Guard: a kept cell type can become near-constant within this subset
  # (e.g., neutrophils nearly absent in all controls), which makes CellDMC
  # rank-deficient. Drop and renormalize.
  v_sub <- apply(frac_sub, 2, var)
  good  <- names(v_sub)[v_sub >= 1e-5]
  if (length(good) < ncol(frac_sub)) {
    cat("  Dropping near-constant-in-subset cell types:",
        paste(setdiff(colnames(frac_sub), good), collapse = ", "), "\n")
  }
  frac_sub <- frac_sub[, good, drop = FALSE]
  frac_sub <- frac_sub / rowSums(frac_sub)
  report_ct_sub <- intersect(report_ct, colnames(frac_sub))
  
  res <- tryCatch(
    CellDMC(beta.m = beta_sub, pheno.v = pheno_v, frac.m = frac_sub,
            adjPMethod = "fdr", adjPThresh = 0.05, mc.cores = 1),
    error = function(e) {
      cat("  CellDMC failed:", conditionMessage(e), "\n"); NULL
    }
  )
  if (is.null(res)) next
  
  # dmct: -1/0/1 per cell type at FDR<0.05
  dmct_counts <- colSums(abs(res$dmct[, -1, drop = FALSE]) > 0, na.rm = TRUE)
  cat("Cell-type-specific DMCs (FDR<0.05):\n"); print(dmct_counts)
  celldmc_summary[[cname]] <- dmct_counts
  
  # Full per-CpG tables for Epi / Neutro, annotated to genes
  for (ct in report_ct_sub) {
    coe <- as.data.frame(res$coe[[ct]]) %>%
      rownames_to_column("probeID") %>%
      left_join(probe_anno_final, by = "probeID") %>%
      arrange(adjP)
    write_csv(coe, file.path(TBL_DIR,
                             paste0("celldmc_", cname, "_", ct, ".csv")))
    cat(sprintf("  %s: %d CpGs at FDR<0.05 (table written)\n",
                ct, sum(coe$adjP < 0.05, na.rm = TRUE)))
  }
}

# Compact summary across contrasts
if (length(celldmc_summary) > 0) {
  cat("\nCellDMC DMC counts (FDR<0.05) across contrasts:\n")
  print(celldmc_summary)
}

# -----------------------------------------------------------------------------
# 15c. CellDMC validation diagnostics — IC vs non-IC contrast
# -----------------------------------------------------------------------------
# Run AFTER section 15b. The IC contrast returned 118 / 235 / 83 DMCs
# (Epi / Neutro / Rest) at FDR<0.05 from only 5 IC cases. Three checks decide
# whether that is real:
#   (1) Covariate-adjusted rerun (age + sex + slide): how many survive?
#   (2) Permutation null: is FDR calibrated at this n, or do shuffled labels
#       also produce tens-to-hundreds of "DMCs"?
#   (3) Leave-one-out over the 5 cases: is one kid carrying the whole signal?
#
# RUNTIME WARNING: the permutation null reruns CellDMC genome-wide B times and
# is the slow part. Start with a small B for a sanity look; raise B and bump
# PERM_CORES for a number you would actually report.

ic_pheno <- contrasts_list[["IC_vs_nonIC"]](ss)
names(ic_pheno) <- rownames(frac_cdmc)          # aligned to beta_collapsed cols
stopifnot(!any(is.na(ic_pheno)))                # IC flag is defined for all 24

# Pull the observed (unadjusted) counts from 16b, recomputing if absent
obs_counts <- celldmc_summary[["IC_vs_nonIC"]]
if (is.null(obs_counts)) {
  res_obs <- CellDMC(beta_collapsed, ic_pheno, frac_cdmc,
                     adjPMethod = "fdr", adjPThresh = 0.05, mc.cores = 1)
  obs_counts <- colSums(abs(res_obs$dmct[, -1, drop = FALSE]) > 0, na.rm = TRUE)
}


## ---- (1) Covariate-adjusted rerun -------------------------------------------
# cov.mod must NOT carry an intercept: frac.m already supplies the per-cell-type
# intercepts, and an all-ones column is collinear with rowSums(frac) = 1, which
# would itself trigger the rank-deficiency error. So we drop column 1.
#
# NOTE: Sex is NA for a few samples (the na_if() in section 2), one of which is
# an IC case, so the adjusted model runs on the complete-covariate subset. To
# keep all 24 instead, drop Sex: model.matrix(~ Age + Slide, ...)[ , -1].
cat("\n--- (1) Covariate-adjusted CellDMC (age + sex + slide) ---\n")

ss$Slide <- factor(substr(ss$Sample_Name, 1, 12))
cov_df   <- data.frame(Age   = suppressWarnings(as.numeric(ss$Age)),
                       Sex   = factor(ss$Sex),
                       Slide = ss$Slide)
cc <- complete.cases(cov_df)
cat("  Complete-covariate samples:", sum(cc), "of", length(cc),
    "| IC cases retained:", sum(ic_pheno[cc] == 1L),
    "| controls:", sum(ic_pheno[cc] == 0L), "\n")

if (sum(ic_pheno[cc] == 1L) >= 4 && sum(ic_pheno[cc] == 0L) >= 4) {
  cov_df_cc <- droplevels(cov_df[cc, ])
  cov_mod   <- model.matrix(~ Age + Sex + Slide, data = cov_df_cc)[, -1, drop = FALSE]
  fr_cc     <- frac_cdmc[cc, , drop = FALSE]; fr_cc <- fr_cc / rowSums(fr_cc)
  
  res_adj <- tryCatch(
    CellDMC(beta.m = beta_collapsed[, cc, drop = FALSE],
            pheno.v = ic_pheno[cc], frac.m = fr_cc, cov.mod = cov_mod,
            adjPMethod = "fdr", adjPThresh = 0.05, mc.cores = 1),
    error = function(e) { cat("  adjusted CellDMC failed:", conditionMessage(e), "\n"); NULL })
  
  if (!is.null(res_adj)) {
    adj_counts <- colSums(abs(res_adj$dmct[, -1, drop = FALSE]) > 0, na.rm = TRUE)
    cat("  Unadjusted vs adjusted DMC counts (FDR<0.05):\n")
    print(rbind(unadjusted = obs_counts[names(adj_counts)],
                adjusted    = adj_counts))
    for (ct in report_ct) {
      coe <- as.data.frame(res_adj$coe[[ct]]) %>%
        rownames_to_column("probeID") %>%
        left_join(probe_anno_final, by = "probeID") %>% arrange(adjP)
      write_csv(coe, file.path(TBL_DIR, paste0("celldmc_IC_adjusted_", ct, ".csv")))
    }
  }
} else {
  cat("  Skipped: too few complete-covariate samples in a group for a stable fit.\n")
}


## ---- (2) Permutation null ---------------------------------------------------
cat("\n--- (2) Permutation null (shuffle IC labels) ---\n")
B_PERM     <- 100      # raise (e.g. 500-1000) for a reportable null
PERM_CORES <- 8        # increase for a faster run
set.seed(20260101)

perm_mat <- matrix(NA_integer_, nrow = B_PERM, ncol = ncol(frac_cdmc),
                   dimnames = list(NULL, colnames(frac_cdmc)))
for (b in seq_len(B_PERM)) {
  perm_pheno <- sample(ic_pheno)               # preserves the 5 / 19 split
  rp <- tryCatch(
    CellDMC(beta_collapsed, perm_pheno, frac_cdmc,
            adjPMethod = "fdr", adjPThresh = 0.05, mc.cores = PERM_CORES),
    error = function(e) NULL)
  if (!is.null(rp))
    perm_mat[b, ] <- colSums(abs(rp$dmct[, -1, drop = FALSE]) > 0, na.rm = TRUE)
  if (b %% 25 == 0) cat("  permutation", b, "/", B_PERM, "\n")
}

perm_ok <- perm_mat[complete.cases(perm_mat), , drop = FALSE]
perm_summary <- map_dfr(colnames(perm_mat), function(ct) {
  pc <- perm_ok[, ct]
  tibble(cell_type   = ct,
         observed    = obs_counts[ct],
         perm_median = median(pc),
         perm_mean   = mean(pc),
         perm_p95    = unname(quantile(pc, 0.95)),
         perm_max    = max(pc),
         # +1 smoothing so the empirical p is never exactly 0
         emp_p_value = (sum(pc >= obs_counts[ct]) + 1) / (length(pc) + 1))
})
cat("  Permutations completed:", nrow(perm_ok), "/", B_PERM, "\n")
print(perm_summary %>% mutate(across(where(is.numeric), ~round(.x, 3))))
write_csv(perm_summary, file.path(TBL_DIR, "celldmc_IC_permutation_null.csv"))
# Read it like this: if perm_median / perm_p95 sit in the tens-to-hundreds and
# emp_p_value is not small, the FDR is NOT calibrated at this n and most of the
# observed hits are noise. If permuted counts are ~0 and emp_p_value is small,
# the signal is above chance.


## ---- (3) Leave-one-out over the IC cases ------------------------------------
cat("\n--- (3) Leave-one-out over the 5 IC cases ---\n")
ic_idx <- which(ic_pheno == 1L)
loo <- map_dfr(ic_idx, function(i) {
  keep <- setdiff(seq_along(ic_pheno), i)
  fr   <- frac_cdmc[keep, , drop = FALSE]; fr <- fr / rowSums(fr)
  rr   <- tryCatch(
    CellDMC(beta_collapsed[, keep, drop = FALSE], ic_pheno[keep], fr,
            adjPMethod = "fdr", adjPThresh = 0.05, mc.cores = 1),
    error = function(e) NULL)
  cc2 <- if (is.null(rr)) setNames(rep(NA_integer_, ncol(frac_cdmc)), colnames(frac_cdmc))
  else colSums(abs(rr$dmct[, -1, drop = FALSE]) > 0, na.rm = TRUE)
  tibble(dropped_sample = names(ic_pheno)[i],
         Epi = cc2["Epi"], Neutro = cc2["Neutro"], Rest = cc2["Rest"])
})
cat("  Full-set counts:",
    paste(names(obs_counts), obs_counts, sep = "=", collapse = ", "), "\n")
print(loo)
write_csv(loo, file.path(TBL_DIR, "celldmc_IC_leave_one_out.csv"))
# If dropping any single case sends a count from ~235 to single/low-double
# digits, that case is driving the result and it is not a stable cohort signal.

# =============================================================================
# 15d-i. Table: CellDMC IC-contrast counts (covariate adjustment + permutation)
# =============================================================================
# Rows are cell types; columns give the
# unadjusted and covariate-adjusted DMC counts and the permutation-null summary.
# "Unadjusted DMCs" is the observed count the permutation p-value tests against.
ct_levels  <- c("Epi", "Neutro", "Rest")
ct_display <- c(Epi = "Epithelial", Neutro = "Neutrophil", Rest = "Rest (pooled)")
adj_vec <- if (exists("adj_counts")) adj_counts[ct_levels] else setNames(rep(NA_integer_, 3), ct_levels)

celldmc_tbl <- perm_summary %>%
  mutate(cell_type = factor(cell_type, levels = ct_levels)) %>%
  arrange(cell_type) %>%
  transmute(
    `Cell type`            = ct_display[as.character(cell_type)],
    `Unadjusted DMCs`      = as.integer(observed),
    `Adjusted DMCs`        = as.integer(adj_vec[as.character(cell_type)]),
    `Null median`          = round(perm_median),
    `Null 95th percentile` = round(perm_p95),
    `Null maximum`         = round(perm_max),
    `Empirical p-value`    = sprintf("%.3f", emp_p_value))

n_perm <- if (exists("perm_ok")) nrow(perm_ok) else NA_integer_
save_pub_table(
  celldmc_tbl,
  file.path(TBL_DIR, "Table_CellDMC_IC_validation"),
  caption = paste0("Table. Cell-type-specific differential methylation for the ",
                   "immunocompromised contrast: covariate adjustment and ",
                   "permutation null."),
  note = paste0("DMC counts are CpGs at FDR < 0.05 from CellDMC. Unadjusted ",
                "DMCs are the base contrast; adjusted DMCs add age, sex, and ",
                "slide covariates. The permutation null shuffled the ",
                "immunocompromised labels ",
                ifelse(is.na(n_perm), "B", n_perm), " times; the empirical ",
                "p-value is the fraction of permutations with a DMC count at ",
                "least as large as the unadjusted observed value. No cell type ",
                "exceeded chance after correction. Rest = pooled non-target ",
                "cells."))


# =============================================================================
# 15d-ii. Figure: CellDMC IC-contrast leave-one-out
# =============================================================================
# Leave-one-out bar plot: DMC count vs the subject dropped, per cell type,
# with a dashed line at the full-cohort count.
loo_long <- loo %>%
  tidyr::pivot_longer(c(Epi, Neutro, Rest),
                      names_to = "cell_type", values_to = "n_dmc") %>%
  mutate(cell_type = dplyr::recode(cell_type, Epi = "Epithelial",
                                   Neutro = "Neutrophil", Rest = "Other"),
         cell_type = factor(cell_type,
                            levels = c("Epithelial", "Neutrophil", "Other")))
full_lines <- tibble(
  cell_type = factor(c("Epithelial", "Neutrophil", "Other"),
                     levels = c("Epithelial", "Neutrophil", "Other")),
  full = c(obs_counts["Epi"], obs_counts["Neutro"], obs_counts["Rest"]))

pC <- ggplot(loo_long, aes(dropped_sample, n_dmc)) +
  geom_col(width = 0.7, fill = "#56B4E9") +
  geom_hline(data = full_lines, aes(yintercept = full),
             linetype = "dashed", colour = "grey30") +
  facet_wrap(~ cell_type, scales = "free_y") +
  labs(x = "Immunocompromised subject removed", y = "DMCs (FDR < 0.05)") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

suppC <- pC +
  labs(title = "CellDMC leave-one-out, immunocompromised contrast",
       caption = paste0("Dashed lines mark full-cohort DMC counts. Dropping a ",
                        "single immunocompromised subject (ARDS_026) inflates ",
                        "epithelial and Rest counts by one to two orders of ",
                        "magnitude, indicating a signal driven by one subject ",
                        "rather than a stable cohort effect. Rest = pooled ",
                        "non-target cells.")) +
  theme(plot.caption = element_text(size = 8, colour = "grey30", hjust = 0))

ggsave(file.path(FIG_DIR, "SuppFig_CellDMC_IC_leave_one_out.pdf"), suppC,
       width = 12, height = 5, device = cairo_pdf, dpi = 600)
ggsave(file.path(FIG_DIR, "SuppFig_CellDMC_IC_leave_one_out.png"), suppC,
       width = 12, height = 5, dpi = 600)
cat("CellDMC IC leave-one-out figure written: SuppFig_CellDMC_IC_leave_one_out.pdf / .png\n")
# -----------------------------------------------------------------------------
# 16. Wrap up + session info
# -----------------------------------------------------------------------------

# =============================================================================
# Save one bundle of the whole analysis (curated; skips objects not in memory)
# =============================================================================
analysis_keep <- c(
  "sample_sheet", "cell_props", "cell_cols", "prop_for_adj",
  "m_collapsed", "beta_collapsed", "m_adjusted", "m_adj_top5",
  "pca_raw", "pca_adj", "var_exp_raw", "var_exp_adj", "variance_compare",
  "cor_table", "sample_hc_adj", "chosen_k", "sil_table",
  "probe_anno_final", "dmc_res", "ora_dir_list",
  "celldmc_summary", "obs_counts", "adj_counts",
  "perm_summary", "perm_ok", "loo")
present <- intersect(analysis_keep, ls(envir = .GlobalEnv))
analysis_bundle <- mget(present, envir = .GlobalEnv)
analysis_bundle$.saved_at <- Sys.time()
analysis_bundle$.session  <- sessionInfo()
saveRDS(analysis_bundle,
        file.path(RDS_DIR, "PARDS_methylome_full_analysis.rds"))
cat("Saved analysis bundle (", length(present), "objects ) to",
    file.path(RDS_DIR, "PARDS_methylome_full_analysis.rds"), "\n")


cat("\n===================================================\n")
cat("Pipeline complete:", format(Sys.time()), "\n")
cat("Figures:  ", FIG_DIR, "\n")
cat("Tables:   ", TBL_DIR, "\n")
cat("RDS data: ", RDS_DIR, "\n")
cat("Log file: ", log_file, "\n")
cat("===================================================\n\n")

print(sessionInfo())
sink()
