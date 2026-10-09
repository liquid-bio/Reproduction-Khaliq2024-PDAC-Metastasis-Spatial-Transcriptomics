#!/bin/bash
# MISTy on all sections (~5 min with 16 cores). Submit after RCTD finishes:
#   sbatch --dependency=afterok:<rctd_job_id> slurm/misty.sbatch
#SBATCH --job-name=misty
#SBATCH --partition=<partition>
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=logs/misty.log
Rscript R/03_misty.R
