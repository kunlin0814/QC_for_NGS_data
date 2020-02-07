#!/bin/bash
#PBS -N CMT-2
#PBS -l walltime=360:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=40gb
#PBS -q batch

cd /scratch/kh31516/Original_Mammary/results/CMT-100/
cat SRR7780979_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
cat SRR7780979_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $2,$3}' >>
mean=$()
> SRR7780979_DepthofCoverage_Distribution.txt
cat SRR7780979_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc