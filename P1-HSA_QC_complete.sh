#!/bin/bash
#PBS -q batch
#PBS -N new-bwa-java-P1
#PBS -l nodes=1:ppn=1
#PBS -l walltime=50:00:00
#PBS -l mem=55gb



# WGS 100 bps 
# WES 76 bps

source_folder='/scratch/kh31516/Pan_cancer/HSA/results/P1' #
script='/work/szlab/kh31516_Lab_Share_script'
file_output='/scratch/kh31516/Pan_cancer/HSA/results/P1'
summary_output='/scratch/kh31516/Pan_cancer/HSA/results/'
DepthOfCoverage='/scratch/kh31516/Pan_cancer/HSA/DepthOfCoverage/P1'
Mutect='/scratch/kh31516/Pan_cancer/HSA/Mutect/P1'
Sequence_length=75
Cancer_Type='HSA'
Normal_Sample='SRR6277192'
Tumor_Sample='SRR6277190'

module load SAMtools/1.9-foss-2016b
module load Anaconda3/2018.12
source activate py35

## Convert the bam file into sam file
cd ${source_folder}/

#samtools view $source_folder/${Normal_Sample}.bam > $file_output/${Normal_Sample}.sam # normal
#samtools view $source_folder/${Tumor_Sample}.bam > $file_output/${Tumor_Sample}.sam # tumor

###### Analyze the BWA mapping and CDS region mapping result ######

cd ${file_output}/

java -cp ${script}/ Line_by_line_Total_dict_Get_exon_reads \
${file_output}/${Normal_Sample}.sam \
${Normal_Sample} \
${file_output}/Normal-${Normal_Sample}-CDS_Mapping_summary.txt \
$Sequence_length

java -cp ${script}/ Line_by_line_Total_dict_Get_exon_reads \
${file_output}/${Tumor_Sample}.sam \
${Tumor_Sample} \
${file_output}/Tumor-${Tumor_Sample}-CDS_Mapping_summary.txt \
$Sequence_length

cat ${file_output}/Normal-${Normal_Sample}-CDS_Mapping_summary.txt >> $summary_output/Total_WES_BWA_CDS_${Cancer_Type}_normal.txt
cat ${file_output}/Tumor-${Tumor_Sample}-CDS_Mapping_summary.txt >>  $summary_output/Total_WES_BWA_CDS_${Cancer_Type}_tumor.txt

##### Sequence Reads Mapping Quality #####

cd ${file_output}/

cat ${file_output}/${Normal_Sample}.sam | cut -f5 > ${file_output}/${Normal_Sample}-mapping_quality

Ngt30=$(cat ${Normal_Sample}-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Ngt60=$(cat ${Normal_Sample}-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ntotal=$(cat ${Normal_Sample}-mapping_quality| wc -l)
Nfra30=$(echo "$((Ngt30))/$((Ntotal))" | bc -l)
Nfra60=$(echo "$((Ngt60))/$((Ntotal))" | bc -l)

printf  "%s\t%4f\t%4f\n" "${Normal_Sample}" "${Nfra30}" "${Nfra60}" >> $summary_output/WES_Total_Mapping_quality_Normal_${Cancer_Type}.txt

cat ${file_output}/${Tumor_Sample}.sam | cut -f5 > ${file_output}/${Tumor_Sample}-mapping_quality

Tgt30=$(cat ${file_output}/${Tumor_Sample}-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Tgt60=$(cat ${file_output}/${Tumor_Sample}-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ttotal=$(cat ${file_output}/${Tumor_Sample}-mapping_quality| wc -l)
Tfra30=$(echo "$((Tgt30))/$((Ttotal))" | bc -l)
Tfra60=$(echo "$((Tgt60))/$((Ttotal))" | bc -l)

printf  "%s\t%4f\t%4f\n" "${Tumor_Sample}" "${Tfra30}" "${Tfra60}" >> $summary_output/WES_Total_Mapping_quality_Tumor_${Cancer_Type}.txt

###### QC of Depth of Coverage and Randomness ######

cd ${DepthOfCoverage}/

cat ${DepthOfCoverage}/${Normal_Sample}_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${DepthOfCoverage}/${Normal_Sample}_DepthofCoverage_Distribution.txt
cat ${DepthOfCoverage}/${Tumor_Sample}_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${DepthOfCoverage}/${Tumor_Sample}_DepthofCoverage_Distribution.txt

#cat ${Tumor_Sample}_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
#cat ${Normal_Sample}_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
python ${script}/Randomness.py ${DepthOfCoverage}/${Normal_Sample}_DepthofCoverage_Distribution.txt ${Normal_Sample} ${Cancer_Type} Normal
python ${script}/Randomness.py ${DepthOfCoverage}/${Tumor_Sample}_DepthofCoverage_Distribution.txt ${Tumor_Sample} ${Cancer_Type} Tumor

cat ${DepthOfCoverage}/${Normal_Sample}_randomness_summary.txt >> $summary_output/Total_WES_${Cancer_Type}_normal_Randomness_Summary.txt
cat ${DepthOfCoverage}/${Tumor_Sample}_randomness_summary.txt >> $summary_output/Total_WES_${Cancer_Type}_tumor_Randomness_Summary.txt

###### Callable bases ########


reference='/work/szlab/Lab_shared_PanCancer/source'
#Yuan_script='/work/szlab/Lab_shared_PanCancer/script'

#cd ${file_output}

#callable=$(cat $Mutect/P1_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1"| wc -l)

#printf  "%s\t%d\n" "i_03A6" "${callable}"  >> $summary_output/Total_WES_Callable_${Cancer_Type}.txt


#rm ${file_output}/${Normal_Sample}.sam
#rm ${file_output}/${Tumor_Sample}.sam
rm ${file_output}/${Normal_Sample}_BWA_summarize
rm ${file_output}/${Tumor_Sample}_BWA_summarize
rm ${file_output}/Normal-${Normal_Sample}-CDS_Mapping_summary.txt
rm ${file_output}/Tumor-${Tumor_Sample}-CDS_Mapping_summary.txt
rm ${file_output}/${Normal_Sample}-mapping_quality ${file_output}/${Tumor_Sample}-mapping_quality
rm ${DepthOfCoverage}/${Normal_Sample}_DepthofCoverage_Distribution.txt
rm ${DepthOfCoverage}/${Tumor_Sample}_DepthofCoverage_Distribution.txt
rm ${DepthOfCoverage}/${Tumor_Sample}_randomness_summary.txt
rm ${DepthOfCoverage}/${Normal_Sample}_randomness_summary.txt