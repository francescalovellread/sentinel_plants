#!/bin/bash

#SBATCH --clusters=all
#SBATCH --partition=short
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --mem-per-cpu=8G 
#SBATCH --job-name=Figs4BCD5ABC
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francesca.lovell-read@merton.ox.ac.uk
#SBATCH --output=out_ex.out
#SBATCH --array=1-50

module load MATLAB #%load MATLAB (can specify version I think if need be)
matlab -nodisplay -nosplash -r "Figs_4BCD_5ABC_RUN($SLURM_ARRAY_TASK_ID); quit;"