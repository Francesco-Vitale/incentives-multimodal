#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --output=ParallelOut

module load scicomp-python-env # use the normal scicomp environment for python
source ~/myvenv/bin/activate
export AMPLPATH=/home/vitalef2/ampl.linux-intel64
srun python master_gurobi.py
