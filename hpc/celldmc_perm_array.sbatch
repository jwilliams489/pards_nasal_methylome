#!/bin/bash
# =============================================================================
# CellDMC permutation null, as a SLURM array job.
#
# 20 tasks x 50 permutations = 1000. Each task regenerates all 1000 permutation
# label vectors from the shared seed and computes only its own slice, so a slice
# is identical no matter which task runs it. Tasks are idempotent: a slice whose
# output already exists is skipped, so resubmitting after a partial failure
# reruns only what is missing.
#
#   sbatch celldmc_perm_array.sbatch          # submit
#   squeue -u "$USER"                         # monitor
#   tail -f "$CELLDMC_WORK/logs"/perm_*_1.out
#
# SITE CONFIGURATION. The three settings below are specific to the cluster this
# was run on and will need changing elsewhere: the partition name, the way the R
# environment is activated, and CELLDMC_WORK. Everything else is portable.
# =============================================================================
#SBATCH --job-name=celldmc_perm
#SBATCH --array=1-20
#SBATCH --cpus-per-task=28
#SBATCH --mem=60G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --exclusive
# #SBATCH --partition=CHANGE_ME       # site-specific; uncomment and set
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=you@example.edu
#
# Add %10 to the array spec (--array=1-20%10) to cap concurrency if the queue
# policy or fairshare requires it.

set -euo pipefail

# Working directory: holds celldmc_perm_input.rds, receives perm_slices/ and
# logs/. Defaults to the directory the job was submitted from.
CELLDMC_WORK="${CELLDMC_WORK:-${SLURM_SUBMIT_DIR:-$PWD}}"
export CELLDMC_WORK
mkdir -p "$CELLDMC_WORK/logs"

# --- Activate the R environment -------------------------------------------
# Replace this block with whatever your site uses. Two common forms:
#
#   module load R/4.5.0
#
#   eval "$(micromamba shell hook --shell bash)" && micromamba activate celldmc
#
if [ -n "${CELLDMC_ENV_INIT:-}" ]; then
  # shellcheck disable=SC1090
  source "$CELLDMC_ENV_INIT"
fi

# Keep BLAS serial; parallelism here is across forked workers, not within them.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$CELLDMC_WORK"

echo "Array task ${SLURM_ARRAY_TASK_ID} on $(hostname), $(date)"
echo "Work dir: $CELLDMC_WORK"
echo "Cores:    ${SLURM_CPUS_PER_TASK}"

Rscript "$(dirname "$0")/celldmc_perm_array.R"

echo "Task ${SLURM_ARRAY_TASK_ID} exited $?, $(date)"
