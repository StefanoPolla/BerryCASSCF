#!/bin/bash
# FULVENE STAGE 2 of 2 -- Berry phase on loops around the intersection located in stage 1.
# Requires a completed stage-1 scan; pass its path with SCAN=...
#SBATCH --job-name=fulvene-berry
#SBATCH --output=logs/fulvene-berry-%j.out
#SBATCH --error=logs/fulvene-berry-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           #### SITE ####
#SBATCH --mem=64G                    #### SITE ####
#SBATCH --time=12:00:00              #### SITE ####
#SBATCH --partition=PARTITION_NAME   #### SITE ####
#SBATCH --account=ACCOUNT_NAME       #### SITE ####

set -euo pipefail

#### SITE #### module loads
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
NPOINTS="${NPOINTS:-17 25}"
RADIUS="${RADIUS:-0.06 15.0}"

# Path to the stage-1 .npz. Override with SCAN=... ; otherwise the newest one is used.
SCAN="${SCAN:-$(ls -t results/fulvene/fulvene_scan_*.npz 2>/dev/null | head -1)}"
if [ -z "${SCAN}" ]; then
    echo "ERROR: no stage-1 scan found. Submit slurm/fulvene_scan.sbatch first," >&2
    echo "       or pass SCAN=path/to/fulvene_scan_*.npz" >&2
    exit 1
fi

echo "host=$(hostname) job=${SLURM_JOB_ID:-none} threads=${OMP_NUM_THREADS}"
echo "scan=${SCAN} basis=${BASIS} cas=${CAS} npoints=${NPOINTS} radius=${RADIUS}"

srun python examples/run_fulvene_berry.py \
    --scan "${SCAN}" \
    --basis "${BASIS}" \
    --cas "${CAS}" \
    --npoints ${NPOINTS} \
    --radius ${RADIUS} \
    --quiet
