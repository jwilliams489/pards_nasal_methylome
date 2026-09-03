# =============================================================================
# celldmc_perm_array.R
#
# One slice of the CellDMC permutation null, run as a SLURM array task.
#
# Each task regenerates ALL B permutation label vectors from the same seed and
# then computes only its own slice. This is what keeps the array reproducible:
# the labels for permutation b are identical regardless of which task computes
# it, or whether it was computed here or by the local loop in
# Methyl_PARDS_pipeline.R (which sets the same RNG kind and seed).
#
# Every analysis parameter -- B, the seed, the FDR threshold, the RNG kind, and
# the observed counts -- travels inside celldmc_perm_input.rds, written by
# PART V of the pipeline. Nothing analytic is set here, so a slice can never be
# computed under different settings than the run it belongs to.
#
# Usage (invoked by celldmc_perm_array.sbatch):
#   Rscript celldmc_perm_array.R
# Reads SLURM_ARRAY_TASK_ID, SLURM_CPUS_PER_TASK and CELLDMC_WORK from the
# environment.
# =============================================================================

suppressPackageStartupMessages({
  library(EpiDISH)
  library(parallel)
})

WORK_DIR  <- Sys.getenv("CELLDMC_WORK", unset = getwd())
INPUT_RDS <- file.path(WORK_DIR, "celldmc_perm_input.rds")
SLICE_DIR <- file.path(WORK_DIR, "perm_slices")
N_TASKS   <- 20L

if (!file.exists(INPUT_RDS))
  stop("Input not found: ", INPUT_RDS,
       "\n  It is written by PART V of Methyl_PARDS_pipeline.R as",
       "\n  analysis_output_v2/celldmc_perm_input.rds. Copy it here.")

dir.create(SLICE_DIR, showWarnings = FALSE, recursive = TRUE)

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
stopifnot(!is.na(task_id), task_id >= 1L, task_id <= N_TASKS)

slice_file <- file.path(SLICE_DIR, sprintf("perm_slice_%03d.rds", task_id))

# Idempotent: a resubmitted array skips slices that already completed.
if (file.exists(slice_file)) {
  cat("Slice", task_id, "already exists; nothing to do.\n")
  quit(save = "no", status = 0)
}

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")

# -----------------------------------------------------------------------------
# Data and parameters, all from the input bundle
# -----------------------------------------------------------------------------
inp <- readRDS(INPUT_RDS)
beta_collapsed <- inp$beta_collapsed
frac_cdmc      <- inp$frac_cdmc
ic_pheno       <- inp$ic_pheno
B_PERM         <- inp$b_perm
SEED           <- inp$seed
FDR_THRESH     <- inp$fdr_thresh
RNG_KIND       <- inp$rng_kind
rm(inp); invisible(gc())

stopifnot(identical(colnames(beta_collapsed), rownames(frac_cdmc)),
          identical(names(ic_pheno), colnames(beta_collapsed)),
          !any(is.na(ic_pheno)),
          is.numeric(B_PERM), is.numeric(SEED), is.numeric(FDR_THRESH))

cat("Task", task_id, "of", N_TASKS, "| cores:", n_cores,
    "| started", format(Sys.time()), "\n")
cat("Node:", Sys.info()[["nodename"]], "\n")
cat("EpiDISH:", as.character(packageVersion("EpiDISH")),
    "| limma:", as.character(packageVersion("limma")), "\n")
cat("B =", B_PERM, "| seed =", SEED, "| FDR <", FDR_THRESH,
    "| RNG =", RNG_KIND, "\n")

dmc_counts <- function(pheno_v) {
  res <- CellDMC(beta.m = beta_collapsed, pheno.v = pheno_v,
                 frac.m = frac_cdmc, adjPMethod = "fdr",
                 adjPThresh = FDR_THRESH, mc.cores = 1)
  colSums(abs(res$dmct[, -1, drop = FALSE]) > 0, na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# Permutation labels: all B generated from one seed, then sliced
# -----------------------------------------------------------------------------
# L'Ecuyer-CMRG gives forked workers independent substreams. The labels are
# drawn in the parent before any forking, so this only guards RNG use inside
# workers -- but the kind must still match what the pipeline's local loop uses,
# or permutation b would differ between the two and the slices could not be
# pooled. Both read the kind from the input bundle for exactly that reason.
RNGkind(RNG_KIND)
set.seed(SEED)
perm_labels <- lapply(seq_len(B_PERM), function(i) sample(ic_pheno))
stopifnot(all(vapply(perm_labels, function(p) sum(p == 1L), integer(1)) ==
                sum(ic_pheno == 1L)))

# Contiguous slices; the final task absorbs any remainder.
per_task    <- ceiling(B_PERM / N_TASKS)
slice_start <- (task_id - 1L) * per_task + 1L
slice_end   <- min(task_id * per_task, B_PERM)
if (slice_start > B_PERM) {
  cat("Nothing to do: slice starts past B.\n")
  quit(save = "no", status = 0)
}
slice_idx <- seq.int(slice_start, slice_end)

cat("Slice covers permutations", slice_start, "to", slice_end,
    "(", length(slice_idx), "total )\n\n")

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
t0 <- Sys.time()
res_list <- mclapply(slice_idx, function(b) {
  tryCatch(dmc_counts(perm_labels[[b]]), error = function(e) NULL)
}, mc.cores = n_cores, mc.preschedule = FALSE)

slice_mat <- matrix(NA_integer_, nrow = length(slice_idx),
                    ncol = ncol(frac_cdmc),
                    dimnames = list(NULL, colnames(frac_cdmc)))
for (j in seq_along(slice_idx)) {
  r <- res_list[[j]]
  if (!is.null(r) && length(r) == ncol(frac_cdmc)) slice_mat[j, ] <- r
}

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
n_ok    <- sum(complete.cases(slice_mat))

saveRDS(list(task_id     = task_id,
             perm_idx    = slice_idx,
             slice_mat   = slice_mat,
             n_ok        = n_ok,
             n_cores     = n_cores,
             node        = Sys.info()[["nodename"]],
             elapsed_min = elapsed,
             seed        = SEED,
             rng_kind    = RNG_KIND,
             b_perm      = B_PERM,
             fdr_thresh  = FDR_THRESH),
        slice_file)

cat(sprintf("\nTask %d done: %d / %d permutations succeeded in %.1f min\n",
            task_id, n_ok, length(slice_idx), elapsed))
cat("Wrote:", slice_file, "\n")
cat("Finished:", format(Sys.time()), "\n")
