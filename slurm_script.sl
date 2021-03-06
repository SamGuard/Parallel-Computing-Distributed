#!/bin/sh
#SBATCH --account=cm30225
#SBATCH --partition=teaching
# Name of job (optional)
#SBATCH --job-name=gj_everyone:)
#SBATCH --nodes=1
#SBATCH --output="output.txt"
#SBATCH --error="error.txt"
mpicc ./Parallel-Computing-Distributed/main.c -o main -lm -O3
mpirun -np 138 ./main 5000 5000 0.0001
mpirun -np 176 ./main 5000 5000 0.0001
