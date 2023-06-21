#!/bin/bash
#SBATCH -J motif_intersection
#SBATCH -o ./output/%j.%x.out
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --export=NONE
#SBATCH --ntasks=5
#SBATCH --time=05:00:00
module load slurm_setup
module load julia/1.8.2

mpiexec -n 5 julia --project="." -- ./scripts/main.jl ../input/cluster/motif_intersection.jl
sacct -j $SLURM_JOB_ID --format=jobid,start,end,CPUTime,Elapsed,ExitCode