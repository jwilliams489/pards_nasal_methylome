# =============================================================================
# celldmc_perm_combine.R
#
# Stitches the array task slices into one permutation matrix and produces the
# summary table. Run after the array completes.
#
# The observed DMC counts are read from celldmc_perm_input.rds, which the
# pipeline writes alongside the data the permutations were computed on. They
# are NOT set here.
#
# That distinction matters. An earlier version of this script hard-coded
#
#     OBS_COUNTS <- c(Epi = 118L, Neutro = 235L, Rest = 83L)
#
# and those counts went stale when the probe set was corrected from 298,646 to
# 298,445, after which the observed counts became 125 / 235 / 88. The null
# distribution was still perfectly valid, but every empirical p-value it
# produced was measured against the wrong observation, and nothing caught it
# because the anchor no longer travelled with the data. Reading the counts from
# the input bundle makes that failure impossible: the anchor and the data are
# the same object.
#
# Usage:
#   Rscript celldmc_perm_combine.R
# Reads CELLDMC_WORK from the environment, defaulting to the working directory.
# =============================================================================

WORK_DIR  <- Sys.getenv("CELLDMC_WORK", unset = getwd())
SLICE_DIR <- file.path(WORK_DIR, "perm_slices")
INPUT_RDS <- file.path(WORK_DIR, "celldmc_perm_input.rds")
OUT_CSV   <- file.path(WORK_DIR, "celldmc_IC_permutation_null_B1000.csv")
OUT_RDS   <- file.path(WORK_DIR, "celldmc_IC_permutation_null_B1000_full.rds")

if (!file.exists(INPUT_RDS))
  stop("Input not found: ", INPUT_RDS,
       "\n  The observed counts and B are read from it; see the header.")

inp        <- readRDS(INPUT_RDS)
OBS_COUNTS <- inp$obs_counts
B_PERM     <- inp$b_perm
SEED       <- inp$seed
RNG_KIND   <- inp$rng_kind
stopifnot(!is.null(OBS_COUNTS), !is.null(B_PERM))

cat("Observed counts (from the input bundle): ",
    paste(sprintf("%s = %d", names(OBS_COUNTS), as.integer(OBS_COUNTS)),
          collapse = ", "), "\n")
cat("B =", B_PERM, "| seed =", SEED, "| RNG =", RNG_KIND, "\n\n")

slice_files <- list.files(SLICE_DIR, pattern = "^perm_slice_\\d+\\.rds$",
                          full.names = TRUE)
cat("Slice files found:", length(slice_files), "\n")
if (length(slice_files) == 0)
  stop("No slices in ", SLICE_DIR, ". Has the array run?")

perm_mat <- matrix(NA_integer_, nrow = B_PERM, ncol = length(OBS_COUNTS),
                   dimnames = list(NULL, names(OBS_COUNTS)))

task_log <- data.frame()
for (f in slice_files) {
  s <- readRDS(f)
  # A slice computed under different settings is not poolable with the rest.
  if (!is.null(s$seed) && !identical(as.numeric(s$seed), as.numeric(SEED)))
    stop("Slice ", basename(f), " used seed ", s$seed, ", expected ", SEED)
  if (!is.null(s$rng_kind) && !identical(s$rng_kind, RNG_KIND))
    stop("Slice ", basename(f), " used RNG kind ", s$rng_kind,
         ", expected ", RNG_KIND,
         "\n  Permutation b is a different relabelling under a different RNG",
         "\n  kind, so these slices cannot be pooled.")
  if (!is.null(s$b_perm) && !identical(as.numeric(s$b_perm), as.numeric(B_PERM)))
    stop("Slice ", basename(f), " used B = ", s$b_perm, ", expected ", B_PERM)
  perm_mat[s$perm_idx, ] <- s$slice_mat
  task_log <- rbind(task_log, data.frame(
    task = s$task_id, n_perm = length(s$perm_idx), n_ok = s$n_ok,
    cores = s$n_cores, node = s$node, minutes = round(s$elapsed_min, 1)))
}
if (nrow(task_log) > 0) {
  task_log <- task_log[order(task_log$task), ]
  cat("\nPer-task summary:\n"); print(task_log, row.names = FALSE)
}

perm_ok <- perm_mat[complete.cases(perm_mat), , drop = FALSE]
cat("\nPermutations completed:", nrow(perm_ok), "of", B_PERM, "\n")
if (nrow(perm_ok) < B_PERM)
  cat("NOTE:", B_PERM - nrow(perm_ok), "permutations missing or failed.",
      "Resubmit the array; completed slices are skipped.\n")

min_p <- 1 / (nrow(perm_ok) + 1)
cat("Smallest achievable empirical p:", signif(min_p, 4), "\n\n")

perm_summary <- do.call(rbind, lapply(colnames(perm_mat), function(ct) {
  pc <- perm_ok[, ct]
  data.frame(
    cell_type   = ct,
    observed    = unname(as.integer(OBS_COUNTS[ct])),
    perm_median = median(pc),
    perm_mean   = round(mean(pc), 1),
    perm_p95    = unname(quantile(pc, 0.95)),
    perm_max    = max(pc),
    # Add-one smoothing, as stated in the manuscript Methods.
    emp_p_value = (sum(pc >= OBS_COUNTS[ct]) + 1) / (length(pc) + 1),
    stringsAsFactors = FALSE
  )
}))

print(perm_summary, row.names = FALSE)
write.csv(perm_summary, OUT_CSV, row.names = FALSE)

saveRDS(list(perm_summary = perm_summary, perm_mat = perm_mat,
             obs_counts = OBS_COUNTS, n_completed = nrow(perm_ok),
             min_p = min_p, task_log = task_log,
             seed = SEED, rng_kind = RNG_KIND,
             session = sessionInfo()),
        OUT_RDS)

cat("\nWrote:", OUT_CSV, "\n      ", OUT_RDS, "\n")
cat("\nCopy the .rds into analysis_output_v2/ and re-run the pipeline;\n")
cat("PART V will load it instead of recomputing the null.\n")
