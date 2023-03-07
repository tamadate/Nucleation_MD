#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time 1:00:00


module load intel
g++ -O3 -o MD.out *cpp -std=c++11 -fopenmp

