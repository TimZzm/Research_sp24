#!/bin/bash

#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --time=6:00:00
#SBATCH --job-name='7e66440'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zz3645@nyu.edu
#SBATCH --output='7e66440.out'

module purge

srun /scratch/work/public/singularity/run-dedalus-3.0.0a0.bash python '7e66440.py' 