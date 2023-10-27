#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -o jobs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -t 1-3000
module load julia/1.6.0

INPUT_ARGS=$(sed -n "$(($SGE_TASK_ID+1))"p input_branch_eelgrass.txt)
julia --project=. -t 1 run_experiment_samples.jl $INPUT_ARGS 30 $JOB_ID data