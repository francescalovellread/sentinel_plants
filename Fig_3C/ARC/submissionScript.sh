#!/bin/bash

#SBATCH --clusters=all
#SBATCH --partition=short
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=8G 
#SBATCH --job-name=Fig3C
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francesca.lovell-read@merton.ox.ac.uk
#SBATCH --output=out_ex.out
#SBATCH --array=1-40

module load MATLAB #%load MATLAB (can specify version I think if need be)
matlab -nodisplay -nosplash -r "Fig_3C_RUN($SLURM_ARRAY_TASK_ID); quit;"