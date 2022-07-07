#!/bin/bash

#SBATCH --job-name=<job-name>
#SBATCH --output=<output-file>
#SBATCH --partition=<partition.
#SBATCH --nodes=1 --cpus-per-task=5 --ntasks=1 --mem=15GB
#SBATCH --array=0-90
#SBATCH --time=8-8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>


module load tools/julia/1.5.2

export JULIA_NUM_THREADS=5


echo
pwd; hostname; date

echo "Starting task $SLURM_ARRAY_TASK_ID"

time julia ./src/main_array-svec.jl $SLURM_ARRAY_TASK_ID 5.0

echo "========================================"
