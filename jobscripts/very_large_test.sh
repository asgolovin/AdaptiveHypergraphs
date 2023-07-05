#!/bin/bash
#SBATCH -J very_large_test
#SBATCH -o ./output/%j.%x.out
#SBATCH --get-user-env
#SBATCH --clusters=cm2_tiny
#SBATCH --partition=cm2_tiny
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --time=02:00:00
module load slurm_setup
module load julia/1.8.2

mpiexec -n 1 julia --project="." -- ./scripts/main.jl ../input/cluster/mpi_very_large_test.jl
sacct -j $SLURM_JOB_ID --format=jobid,start,end,CPUTime,Elapsed,ExitCode