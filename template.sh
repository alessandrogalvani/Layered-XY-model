#!/bin/bash
#SBATCH --job-name=clocktutto    # Job name
#SBATCH -N1
#SBATCH --ntasks=20                  
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem-per-cpu=1gb                     # Job memory request
#SBATCH --time=09:00:00               # Time limit hrs:min:sec
#SBATCH --output=out_clockstrati_%j.txt
#SBATCH --partition=regular1,regular2
pwd; hostname; date

#module load gnu/4.9.2
#module load intel/14.0 
#module load intel/18.0.3.222

#g++ clockStrati.cpp
#./a.out
mpiCC -std=c++11 -mcmodel=medium clockStrati.cpp -o Out.out
mpirun Out.out     
