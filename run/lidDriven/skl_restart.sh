#!/bin/bash
#SBATCH --job-name=qE_30

#SBATCH --partition=general
#SBATCH --mem=177000
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time 47:55:00

###SBATCH --mail-type=ALL
###SBATCH --mail-user=daniel.suarez.cambra@upc.edu
  
### Launch the parallel job on the allocated compute nodes
printf "Running Q2DMHDFoam  ...\n"
Q2DmhdFoam_32 > logs/log.q2dmhdFoam
printf "Running foamToVTK ...\n"
foamToVTK -latestTime > logs/log.foamToVtK
