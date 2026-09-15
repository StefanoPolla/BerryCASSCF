#!/bin/bash
#SBATCH --job-name=berry-loop
#SBATCH --output=logs/berry-loop-%j.out
#SBATCH --error=logs/berry-loop-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8            #### SITE #### match to the node's core count
#SBATCH --mem=32G                    #### SITE #### see docs/compute.md for per-CAS estimates
#SBATCH --time=08:00:00              #### SITE ####
#SBATCH --partition=PARTITION_NAME   #### SITE ####
#SBATCH --account=ACCOUNT_NAME       #### SITE ####

set -euo pipefail

#### SITE #### module loads, e.g.:
# module purge
# module load 2023 Python/3.12.3-GCCcore-13.3.0

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$PROJECT_DIR"
mkdir -p logs results/berry

#### SITE #### environment activation
source .venv/bin/activate

# PySCF threading: keep BLAS threads equal to the allocation, no oversubscription.
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${OMP_NUM_THREADS}"
export OPENBLAS_NUM_THREADS="${OMP_NUM_THREADS}"
export PYSCF_MAX_MEMORY=$((${SLURM_MEM_PER_NODE:-32000} * 8 / 10))

# ---------------------------------------------------------------------------
# Arguments. Override on submission:
#   sbatch --export=ALL,LOOPS="C_x",CAS="6,6",NPOINTS="25 41" slurm/berry_loop.sbatch
# ---------------------------------------------------------------------------
LOOPS="${LOOPS:-C_x C_1 C_2}"
CAS="${CAS:-6,6}"
NPOINTS="${NPOINTS:-25 41}"
BASIS="${BASIS:-sto-3g}"

echo "host=$(hostname) job=${SLURM_JOB_ID:-none} threads=${OMP_NUM_THREADS}"
echo "loops=${LOOPS} cas=${CAS} npoints=${NPOINTS} basis=${BASIS}"

srun python examples/run_formaldimine_berry.py \
    --loops ${LOOPS} \
    --cas ${CAS} \
    --npoints ${NPOINTS} \
    --basis "${BASIS}" \
    --quiet

python examples/summarize.py
