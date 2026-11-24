#!/bin/bash
#SBATCH --job-name=pmwDPR
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=pmwDPR_%j.log
set -euo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh" # CHANGE to your conda install
conda activate pmwDPR
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
pmwDPR run examples/gmi_dpr.yaml
