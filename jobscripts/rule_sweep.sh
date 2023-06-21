#!/bin/bash
#SBATCH -J prop_voting_rts_p_sweep
#SBATCH -o ./output/%j.%x.out
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --export=NONE
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=28
#SBATCH --time=20:00:00
module load slurm_setup
module load julia/1.8.2

mpiexec -n 112 julia --project="." -- ./scripts/main.jl ../input/cluster/prop_voting_rts.jl
sacct -j $SLURM_JOB_ID --format=jobid,start,end,CPUTime,Elapsed,ExitCode,MaxRSS