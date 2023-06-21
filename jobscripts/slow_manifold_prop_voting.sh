#!/bin/bash
#SBATCH -J slow_man_prop_rts
#SBATCH -o ./output/%j.%x.out
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --export=NONE
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --time=04:00:00
module load slurm_setup
module load julia/1.8.2

mpiexec -n 56 julia --project="." -- ./scripts/main.jl ../input/cluster/slow_manifold_prop_voting_rts.jl
sacct -j $SLURM_JOB_ID --format=jobid,start,end,CPUTime,Elapsed,ExitCode,MaxRSS