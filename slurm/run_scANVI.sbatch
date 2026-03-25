#!/bin/bash
#SBATCH --job-name=scANVI_training
#SBATCH --output=scANVI_%j.out
#SBATCH --error=scANVI_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

# Use absolute log paths if you like:
# #SBATCH --output=/path/to/integration_pipeline/logs/scANVI_%j.out

export ATLAS_PIPELINE_ROOT="${ATLAS_PIPELINE_ROOT:-/path/to/integration_pipeline}"

source ~/.bashrc
conda activate your_scvi_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd "$ATLAS_PIPELINE_ROOT" || exit 1
python scVI_h_scANVI.py
