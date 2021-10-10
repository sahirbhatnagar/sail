#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -o log/
#PBS -e log/
#PBS -N WeakFit
#PBS -m ea
#PBS -l mem=8G
#PBS -l vmem=8G
#PBS -l nodes=1:ppn=1
#PBS -t 1-200

cd $PBS_O_WORKDIR

SRC=$PBS_O_WORKDIR

Rscript ${SRC}/mcgill_simulations_gendata.R ${index}


#Rscript /mnt/GREENWOOD_BACKUP/home/sahir.bhatnagar/sail/sail_git_v2/sail/my_sims/mcgill_simulations_gendata.R ${index}
