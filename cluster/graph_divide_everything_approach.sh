#!/bin/sh

<<<<<<< HEAD
#SBATCH -a 1-6
=======
#SBATCH -a 1-7
>>>>>>> e7e5eb9683cecb558052fbb4181a18066ec06eb1
#SBATCH --cpus-per-task=2
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH -p sched_mit_sloan_batch
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sgilmour@mit.edu
#SBATCH --output=../logs/graph_divide_everything_approach_\%a.log

module load julia/1.1.0
module load gurobi/8.0.1

srun julia ../scripts/graph_approach/graph_approach_divide_everything.jl $SLURM_ARRAY_TASK_ID
