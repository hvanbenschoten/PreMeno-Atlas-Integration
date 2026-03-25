#!/bin/bash
#SBATCH --job-name=preintegration
#SBATCH --output=preintegration_%j.out
#SBATCH --error=preintegration_%j.err
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --partition=standard

# Set to your clone of this repository (folder containing pipeline_paths.py)
export ATLAS_PIPELINE_ROOT="${ATLAS_PIPELINE_ROOT:-/path/to/integration_pipeline}"

source ~/.bashrc
conda activate your_scanpy_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd "$ATLAS_PIPELINE_ROOT" || exit 1
python preintegration.py
