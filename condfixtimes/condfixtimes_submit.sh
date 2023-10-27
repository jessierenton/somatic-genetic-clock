#!/bin/bash
#$ -pe smp 20
#$ -l h_vmem=7G
#$ -l h_rt=240:0:0
#$ -cwd
#$ -j y
#$ -N condfixtimes

module load julia/1.6.0

julia --project=.. --threads=20 condfixtimes_r8.jl
