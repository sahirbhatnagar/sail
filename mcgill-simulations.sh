#!/bin/bash
#PBS -l walltime=0:40:00
#PBS -o log/
#PBS -e log/
#PBS -N sim2
#PBS -m ea
#PBS -l mem=12G
#PBS -l vmem=12G
#PBS -l nodes=1:ppn=10
#PBS -t 1-500

cd $PBS_O_WORKDIR

SRC=$PBS_O_WORKDIR

Rscript ${SRC}/mcgill_simulations_gendata.R
