#!/bin/bash

#SBATCH --nodes=12
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --time=50:00:00
#SBATCH --job-name='3d_A12_rame6_43'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zz3645@nyu.edu
#SBATCH --output='3d_A12_RaMe6_43.out'

module purge

srun /scratch/work/public/singularity/run-dedalus-3.0.0a0.bash python 'A12_RaMe6_43.py' 