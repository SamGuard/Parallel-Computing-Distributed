#!/bin/sh
#SBATCH --account=cm30225
#SBATCH --partition=teaching
# Name of job (optional)
#SBATCH --job-name=aAAaaAAaa
#SBATCH --nodes=4
#SBATCH --output="output.txt"
#SBATCH --error="error.txt"
mpicc ./Parallel-Computing-Distributed/main.c -o main -lm -O3
mpirun -np 2 ./main 1000 1000 0.0001
mpirun -np 3 ./main 1000 1000 0.0001
mpirun -np 4 ./main 1000 1000 0.0001
mpirun -np 8 ./main 1000 1000 0.0001
mpirun -np 16 ./main 1000 1000 0.0001
