#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=3G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -o jobs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -t 1-90
module load julia/1.6.0

INPUT_ARGS=$(sed -n "$(($SGE_TASK_ID+1))"p input.txt)
julia --project=. -t 1 run_moran_fulldist.jl $INPUT_ARGS 10 0.1 20 split data_tenyrs
