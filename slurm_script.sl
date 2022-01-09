#!/bin/sh
#SBATCH --account=cm30225
#SBATCH --partition=teaching
# Name of job (optional)
#SBATCH --job-name=aAAaaAAaa
#SBATCH --nodes=4
#SBATCH --output="output.txt"
#SBATCH --error="error.txt"
mpicc ./Parallel-Computing-Distributed/main.c -o main -lm -O3
mpirun -np 104 ./main 40000 40000 0.0001
mpirun -np 138 ./main 40000 40000 0.0001
mpirun -np 176 ./main 40000 40000 0.0001
