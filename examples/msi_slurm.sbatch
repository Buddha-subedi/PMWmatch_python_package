#!/bin/bash
#SBATCH --job-name=pmwmatch
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=pmwmatch_%j.log
set -euo pipefail
source "$HOME/miniconda3/etc/profile.d/conda.sh" # CHANGE to your conda install
conda activate pmwmatch
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
pmwmatch run examples/gmi_dpr.yaml
