#!/bin/bash
#PBS -N CMT-2
#PBS -l walltime=350:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=50gb
#PBS -q batch


module load Anaconda3/2019.03
source activate py35