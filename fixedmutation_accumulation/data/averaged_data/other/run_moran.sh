#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -o jobs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -t 1-4
module load julia/1.6.0

ENV=/data/home/hfx987/multilevelselection/molecularclock_maintxt_010323
INPUT_ARGS=$(sed -n "$(($SGE_TASK_ID+1))"p input.txt)
julia --project=$ENV -t 1 run_moran.jl $INPUT_ARGS 50 0.5 20 data