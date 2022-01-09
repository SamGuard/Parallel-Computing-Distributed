#!/bin/sh
#SBATCH --account=cm30225
#SBATCH --partition=teaching
# Name of job (optional)
#SBATCH --job-name=low_workers
#SBATCH --nodes=1
#SBATCH --output="low_output.txt"
#SBATCH --error="low_error.txt"
mpicc ./Parallel-Computing-Distributed/main.c -o main -lm -O3
mpirun -np 3 ./main 5000 5000 0.0001
mpirun -np 4 ./main 5000 5000 0.0001
mpirun -np 8 ./main 5000 5000 0.0001
mpirun -np 16 ./main 5000 5000 0.0001
