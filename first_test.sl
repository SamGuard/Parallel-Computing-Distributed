#!/bin/sh
#SBATCH --account=cm30225
#SBATCH --partition=teaching
# Name of job (optional)
#SBATCH --job-name=dist_test_1
#SBATCH --nodes=1
#SBATCH --output="output.txt"
#SBATCH --error="error.txt"
mpicc ./Parallel-Computing-Distributed/main.c -o main -lm
mpirun -np 44 ./main