#!/bin/bash

#SBATCH --nodes=2
#SBATCH --tasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --time=6:00:00
#SBATCH --job-name='1e76440'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zz3645@nyu.edu
#SBATCH --output='1e76440.out'

module purge

srun /scratch/work/public/singularity/run-dedalus-3.0.0a0.bash python '1e76440.py' 