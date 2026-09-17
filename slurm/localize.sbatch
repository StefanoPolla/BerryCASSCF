#!/bin/bash
#SBATCH --job-name=localize
#SBATCH --output=logs/localize-%j.out
#SBATCH --error=logs/localize-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8            #### SITE #### match to the node's core count
#SBATCH --mem=64G                    #### SITE #### CAS(12,12) butadiene is 853k determinants
#SBATCH --time=48:00:00              #### SITE #### see the estimate below; be generous
#SBATCH --partition=PARTITION_NAME   #### SITE ####
#SBATCH --account=ACCOUNT_NAME       #### SITE ####

set -euo pipefail

#### SITE #### module loads, e.g.:
# module purge
# module load 2023 Python/3.12.3-GCCcore-13.3.0

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$PROJECT_DIR"
mkdir -p logs results/localize

#### SITE #### environment activation
source .venv/bin/activate

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${OMP_NUM_THREADS}"
export OPENBLAS_NUM_THREADS="${OMP_NUM_THREADS}"
export PYSCF_MAX_MEMORY=$((${SLURM_MEM_PER_NODE:-64000} * 8 / 10))

# ---------------------------------------------------------------------------
# Locate what loop transport encircles, by bisection and triangulation.
# docs/todo.md section 9; method explained in notebooks/locating_intersections.ipynb.
#
#   sbatch --export=ALL,SYSTEM=butadiene,CAS="12 12" slurm/localize.sbatch
#
# The driver skips a run whose result file already exists; pass FORCE=--force to redo.
# ---------------------------------------------------------------------------
SYSTEM="${SYSTEM:-butadiene}"
CAS="${CAS:-12 12}"
CENTRES="${CENTRES:-3}"
FORCE="${FORCE:-}"

echo "host=$(hostname) job=${SLURM_JOB_ID:-none} threads=${OMP_NUM_THREADS}"
echo "system=${SYSTEM} cas=${CAS} centres=${CENTRES}"

srun python examples/run_localization.py "${SYSTEM}" \
    --cas ${CAS} \
    --centres ${CENTRES} \
    ${FORCE}
