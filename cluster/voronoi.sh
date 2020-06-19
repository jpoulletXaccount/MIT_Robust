#!/bin/sh

#SBATCH -a 1-7
#SBATCH --cpus-per-task=2
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH -p sched_mit_sloan_batch
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sgilmour@mit.edu
#SBATCH --output=../logs/voronoi_\%a.log

module load julia/1.1.0
module load gurobi/8.0.1

srun julia ../scripts/mip_count_partitions_variables/voronoi.jl $SLURM_ARRAY_TASK_ID
