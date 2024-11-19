#!/bin/bash

#SBATCH --nodes=6
#SBATCH --tasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=128GB
#SBATCH --time=100:00:00
#SBATCH --job-name='3d_A12_n1o24'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zz3645@nyu.edu
#SBATCH --output='3d_A12_n1o24.out'

module purge

srun /scratch/work/public/singularity/run-dedalus-3.0.0a0.bash python 'A12_n1o24.py' 