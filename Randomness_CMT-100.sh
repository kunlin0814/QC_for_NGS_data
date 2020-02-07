#!/bin/bash
#PBS -N CMT-2
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -q batch

cd /scratch/kh31516/Original_Mammary/results/CMT-100/
module load Anaconda3/2019.03
source activate py35

cat SRR7780976_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > SRR7780976-DepthofCoverage_Distribution.txt
#cat SRR7780979_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
python /scratch/kh31516/Original_Mammary/Randomness.py SRR7780976-DepthofCoverage_Distribution.txt SRR7780976

mv SRR7780976_randomness_summary.txt /scratch/kh31516/Original_Mammary/results/Randomness/Normal/
cat /scratch/kh31516/Original_Mammary/results/Randomness/Normal/SRR7780976_randomness_summary.txt >> /scratch/kh31516/Original_Mammary/Total_Mammary_Normal_Randomness_Summary.txt

cat SRR7780979_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > SRR7780979_DepthofCoverage_Distribution.txt
python /scratch/kh31516/Original_Mammary/Randomness.py SRR7780979-DepthofCoverage_Distribution.txt SRR7780979
mv SRR7780979_randomness_summary.txt /scratch/kh31516/Original_Mammary/results/Randomness/Tumor/
cat /scratch/kh31516/Original_Mammary/results/Randomness/Tumor/SRR7780979_randomness_summary.txt >> /scratch/kh31516/Original_Mammary/Total_Mammary_Tumor_Randomness_Summary.txt
