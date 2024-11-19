#!/bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=04:00:00
#SBATCH --job-name='post_1e76440'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zz3645@nyu.edu
#SBATCH --output='postanalysis_1e76440.out'

module purge

srun /scratch/work/public/singularity/run-dedalus-3.0.0a0.bash python 'postanalysis.py' 