#!/bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=04:00:00
#SBATCH --job-name='post_wzc'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zz3645@nyu.edu
#SBATCH --output='wzc_DDP.out'

module purge

srun /scratch/work/public/singularity/run-dedalus-3.0.0a0.bash python 'for_wzc.py' 