#!/bin/bash
#PBS -l walltime=1:00:00
#PBS -o log/
#PBS -e log/
#PBS -N glintLassoBT
#PBS -m ea
#PBS -l mem=8G
#PBS -l vmem=8G
#PBS -l nodes=1:ppn=1


cd $PBS_O_WORKDIR

SRC=$PBS_O_WORKDIR

Rscript ${SRC}/mcgill_simulations_gendata_glint_lassoBT.R ${index}
