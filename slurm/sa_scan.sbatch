#!/bin/bash
#SBATCH --job-name=sa-scan
#SBATCH --output=logs/sa-scan-%j.out
#SBATCH --error=logs/sa-scan-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16           #### SITE ####
#SBATCH --mem=64G                    #### SITE ####
#SBATCH --time=24:00:00              #### SITE ####
#SBATCH --partition=PARTITION_NAME   #### SITE ####
#SBATCH --account=ACCOUNT_NAME       #### SITE ####

set -euo pipefail

#### SITE #### module loads
# module purge
# module load 2023 Python/3.12.3-GCCcore-13.3.0

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$PROJECT_DIR"
mkdir -p logs results/scan

#### SITE #### environment activation
source .venv/bin/activate

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${OMP_NUM_THREADS}"
export OPENBLAS_NUM_THREADS="${OMP_NUM_THREADS}"

# The scan checkpoints after every grid row, so re-submitting this exact job resumes
# rather than restarting. Safe to run under a short time limit and requeue.
LOOPS="${LOOPS:-C_x}"
CAS="${CAS:-6,6}"
GRID="${GRID:-41 41}"
BASIS="${BASIS:-sto-3g}"

echo "host=$(hostname) job=${SLURM_JOB_ID:-none} threads=${OMP_NUM_THREADS}"
echo "loops=${LOOPS} cas=${CAS} grid=${GRID} basis=${BASIS}"

srun python examples/run_formaldimine_scan.py \
    --loops ${LOOPS} \
    --cas ${CAS} \
    --grid ${GRID} \
    --basis "${BASIS}" \
    --quiet
