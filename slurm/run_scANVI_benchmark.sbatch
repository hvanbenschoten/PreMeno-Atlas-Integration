#!/bin/bash
#SBATCH --job-name=scANVI_benchmark
#SBATCH --output=scANVI_benchmark_%j.out
#SBATCH --error=scANVI_benchmark_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

export ATLAS_PIPELINE_ROOT="${ATLAS_PIPELINE_ROOT:-/path/to/integration_pipeline}"

source ~/.bashrc
conda activate your_scvi_env

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd "$ATLAS_PIPELINE_ROOT" || exit 1
python scANVI_benchmark.py
