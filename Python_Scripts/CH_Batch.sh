#!/bin/bash
#SBATCH --job-name=CH_Python
#SBATCH --time=512:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=long

# send an email when each individual array task starts/ends/fails
#SBATCH --mail-user=dmasse@umass.edu
#SBATCH --mail-type=FAIL,ARRAY_TASKS

# Define what values the variable "$SLURM_ARRAY_TASK_ID" should take on:
#SBATCH --array=1-8%200
# The placeholder "%a" will take on the current array value
#SBATCH --output=/home/masse/R_sim_p%a.out

# Run our R script for each value of $SLURM_ARRAY_TASK_ID defined in the ``--array`` argument
# Activate conda environment:
module load miniconda

eval "$(conda shell.bash hook)"

conda activate py_ch

python CH_Batch.py $SLURM_ARRAY_TASK_ID