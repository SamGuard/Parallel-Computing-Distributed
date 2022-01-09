#!/bin/sh
#SBATCH --account=cm30225
#SBATCH --partition=teaching
# Name of job (optional)
#SBATCH --job-name=seq
#SBATCH --nodes=1
#SBATCH --output="seq_output.txt"
#SBATCH --error="seq_error.txt"
gcc ./Parallel-Computing-Distributed/seq_main.c -o seq -O3
./seq 1000 1000 0.0001
./seq 5000 5000 0.0001
./seq 10000 10000 0.0001
./seq 20000 20000 0.0001
./seq 40000 40000 0.0001
