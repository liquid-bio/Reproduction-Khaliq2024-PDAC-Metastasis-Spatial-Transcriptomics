#!/bin/bash
# One RCTD task per section (line N+1 of data/samples.csv). Peak memory ~25 GB, 25-40 min with 16 cores.
# Submit from the repository root: sbatch slurm/rctd_array.sbatch
#SBATCH --job-name=rctd
#SBATCH --partition=<partition>
#SBATCH --array=1-29
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --output=logs/rctd_%a.log
s=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" data/samples.csv | cut -d, -f1)
Rscript R/02_rctd.R "$s"
