#!/bin/bash
#PBS -N CMT-2
#PBS -l walltime=360:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=40gb
#PBS -q batch

cat SRR7780979_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | tail -n +3  awk 'BEGIN{OFS="\t"} 
> SRR7780979_DepthofCoverage_Distribution.txt
cat SRR7780979_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc