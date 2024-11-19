#!/bin/bash

#SBATCH --nodes=2
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=160GB
#SBATCH --time=4:00:00
#SBATCH --job-name='cjt'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zz3645@nyu.edu
#SBATCH --output='cjt.out'

module purge

srun /scratch/work/public/singularity/run-dedalus-3.0.0a0.bash python '1024_256.py' 