# -----------------------------------------------------------------------------
# Re-anchor a saved CellDMC permutation null to the current observed counts.
# -----------------------------------------------------------------------------
# The permutation null is expensive (1,000 genome-wide CellDMC fits), so it is
# computed once on a cluster and reloaded by Methyl_PARDS_pipeline.R. If the
# observed DMC counts later change -- as they did here when the probe set was
# corrected from 298,646 to 298,445 -- the stored empirical p-values are still
# measured against the OLD observed counts and are wrong, even though the null
# distribution itself is perfectly valid.
#
# The null does not need recomputing in that situation. The empirical p-value is
# a function of the permutation counts and the observed count only, so it can be
# recomputed exactly from the stored matrix in a second.
#
# Usage:
#   Rscript reanchor_permutation_null.R <in.rds> <out.rds> <Epi> <Neutro> <Rest>
# Example, for the counts reported in this manuscript:
#   Rscript reanchor_permutation_null.R \
#     celldmc_IC_permutation_null_B1000_full.rds \
#     celldmc_IC_permutation_null_B1000_full_reanchored.rds 125 235 88
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5)
  stop("Usage: reanchor_permutation_null.R <in.rds> <out.rds> <Epi> <Neutro> <Rest>")

in_path  <- args[1]
out_path <- args[2]
new_obs  <- setNames(as.integer(args[3:5]), c("Epi", "Neutro", "Rest"))

x <- readRDS(in_path)
stopifnot(all(c("perm_mat", "perm_summary") %in% names(x)))

pm <- x$perm_mat
pm <- pm[complete.cases(pm), , drop = FALSE]
B  <- nrow(pm)
cat("Permutations available:", B, "\n")
if (B < nrow(x$perm_mat))
  warning("Only ", B, " of ", nrow(x$perm_mat), " permutations completed.")

cell_types <- x$perm_summary$cell_type
stopifnot(all(cell_types %in% colnames(pm)),
          all(cell_types %in% names(new_obs)))

# Add-one smoothing, matching the manuscript Methods: the empirical p-value is
# the proportion of permutations whose count equalled or exceeded the observed
# count, computed as (1 + k) / (1 + B).
emp_p <- vapply(cell_types, function(ct) {
  (1 + sum(pm[, ct] >= new_obs[[ct]])) / (1 + B)
}, numeric(1))

old <- x$perm_summary
cat("\n")
cat(sprintf("%-10s %10s %10s %10s %10s\n",
            "cell_type", "old_obs", "new_obs", "old_p", "new_p"))
for (i in seq_along(cell_types))
  cat(sprintf("%-10s %10d %10d %10.4f %10.4f\n",
              cell_types[i], as.integer(old$observed[i]),
              new_obs[[cell_types[i]]], old$emp_p_value[i], emp_p[i]))

x$perm_summary$observed    <- as.integer(new_obs[cell_types])
x$perm_summary$emp_p_value <- unname(emp_p)
if (!is.null(x$obs_counts)) x$obs_counts <- new_obs[names(x$obs_counts)]

# The null distribution itself is untouched: median, mean, p95 and max describe
# the permutations, not the observation, so they are carried over unchanged.

saveRDS(x, out_path)
cat("\nWrote", out_path, "\n")
cat("Place it in analysis_output_v2/ as celldmc_IC_permutation_null_B1000_full.rds\n")
cat("so the pipeline loads the re-anchored null on the next run.\n")
