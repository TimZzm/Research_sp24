#!/bin/bash

#SBATCH --nodes=2
#SBATCH --tasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --time=4:00:00
#SBATCH --job-name='draft_3d_A4'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zz3645@nyu.edu
#SBATCH --output='draft_3d_A4.out'

module purge

srun /scratch/work/public/singularity/run-dedalus-3.0.0a0.bash python 'Draft_one.py' 