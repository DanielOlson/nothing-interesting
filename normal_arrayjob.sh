#!/bin/bash
#SBATCH -J trf_big # A single job name for the array
#SBATCH -N 1 # All cores on one Node
#SBATCH -p wheeler_lab_large_cpu 
#SBATCH --mem 64G # Memory request
#SBATCH -t 1-2:00 # Maximum execution time (D-HH:MM)

cd /home/place/UPAPE
(/usr/bin/time -v ./trf ./split/t2t_"${SLURM_ARRAY_TASK_ID}".fa 2 7 7 80 10 30 2000 -h) &> ./trf_2000_time/time_"${SLURM_ARRAY_TASK_ID}".txt
321

#usage is 'sbatch --array 0-<max> normal_arrayjob.sh'