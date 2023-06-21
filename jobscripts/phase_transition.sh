#!/bin/bash
#SBATCH -J phase_transition_sweep
#SBATCH -o ./output/%j.%x.out
#SBATCH --get-user-env
#SBATCH --clusters=inter
#SBATCH --partition=cm2_inter
#SBATCH --export=NONE
#SBATCH --nodes=4
#SBATCH --time=00:07:00
module load slurm_setup
module load julia/1.8.5

echo $LD_LIBRARY_PATH
mpiexec -n 96 julia --project="." -- ./scripts/main.jl ../input/cluster/phase_transition_p_small.jl
sacct -j $SLURM_JOB_ID --format=jobid,start,end,CPUTime,Elapsed,ExitCode