#!/bin/bash
#SBATCH -J motifs_all_rules
#SBATCH -o ./output/%j.%x.out
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --export=NONE
#SBATCH --ntasks=24
#SBATCH --time=01:20:00
module load slurm_setup
module load julia/1.8.2

mpiexec -n 24 julia --project="." -- ./scripts/main.jl ../input/cluster/motifs_all_rules.jl
sacct -j $SLURM_JOB_ID --format=jobid,start,end,CPUTime,Elapsed,ExitCode,MaxRSS