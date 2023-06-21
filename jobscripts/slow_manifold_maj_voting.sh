#!/bin/bash
#SBATCH -J slow_man_maj_rtr
#SBATCH -o ./output/%j.%x.out
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=00:40:00
module load slurm_setup
module load julia/1.8.2

mpiexec -n $SLURM_NTASKS julia --project="." -- ./scripts/main.jl ../input/cluster/slow_manifold_maj_voting_rtr.jl
sacct -j $SLURM_JOB_ID --format=jobid,start,end,CPUTime,Elapsed,ExitCode,MaxRSS