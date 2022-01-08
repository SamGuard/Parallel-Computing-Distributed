#!/bin/sh
#SBATCH --account=cm30225
#SBATCH --partition=teaching
# Name of job (optional)
#SBATCH --job-name=aAAaaAAaa
#SBATCH --nodes=4
#SBATCH --output="output.txt"
#SBATCH --error="error.txt"
mpicc ./Parallel-Computing-Distributed/main.c -o main -lm
mpirun ./main 1024 1024 0.001
