#!/bin/bash
# FULVENE STAGE 1 of 2 -- locate the S1/S0 intersection in the (r, theta) plane.
# Run this first. Stage 2 (slurm/fulvene_berry.sbatch) reads its output.
#SBATCH --job-name=fulvene-scan
#SBATCH --output=logs/fulvene-scan-%j.out
#SBATCH --error=logs/fulvene-scan-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           #### SITE ####
#SBATCH --mem=64G                    #### SITE ####
#SBATCH --time=12:00:00              #### SITE #### see docs/compute.md for the estimate
#SBATCH --partition=PARTITION_NAME   #### SITE ####
#SBATCH --account=ACCOUNT_NAME       #### SITE ####

set -euo pipefail

#### SITE #### module loads, e.g.:
# module purge
# module load 2023 Python/3.12.3-GCCcore-13.3.0

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$PROJECT_DIR"
mkdir -p logs results/fulvene

#### SITE #### environment activation
source .venv/bin/activate

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${OMP_NUM_THREADS}"
export OPENBLAS_NUM_THREADS="${OMP_NUM_THREADS}"

BASIS="${BASIS:-6-31g*}"
CAS="${CAS:-6,6}"
GRID="${GRID:-21 21}"

echo "host=$(hostname) job=${SLURM_JOB_ID:-none} threads=${OMP_NUM_THREADS}"
echo "basis=${BASIS} cas=${CAS} grid=${GRID}"

# Checkpoints after every grid row: re-submitting this identical job resumes rather than
# restarting, so a short time limit with requeue is safe.
srun python examples/run_fulvene_scan.py \
    --basis "${BASIS}" \
    --cas "${CAS}" \
    --grid ${GRID} \
    --quiet
