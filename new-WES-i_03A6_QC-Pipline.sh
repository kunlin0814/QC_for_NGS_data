#!/bin/bash
#PBS -q batch
#PBS -N new-bwa-java
#PBS -l nodes=1:ppn=1
#PBS -l walltime=50:00:00
#PBS -l mem=55gb
#PBS -M kh31516@uga.edu 
#PBS -m ae


# WGS 100 bps 
# WES 76 bps

source_folder='/scratch/kh31516/Pan_cancer/glioma/results/WES/i_03A6' #
script='/work/szlab/kh31516_Lab_Share_script'
file_output='/scratch/kh31516/Pan_cancer/glioma/results'
summary_output='/scratch/kh31516/Pan_cancer/glioma/results'
Sequence_length=75

module load SAMtools/1.9-foss-2016b
module load Anaconda3/2018.12
source activate py35

## Convert the bam file into sam file
cd ${source_folder}/

samtools view $source_folder/SRR10351810.bam > $file_output/SRR10351810.sam # normal
samtools view $source_folder/SRR10351811.bam > $file_output/SRR10351811.sam # tumor

###### Analyze the BWA mapping and CDS region mapping result ######

cd ${file_output}/

java -cp ${script}/ Line_by_line_Total_dict_Get_exon_reads \
${file_output}/SRR10351810.sam \
SRR10351810 \
${file_output}/Normal-SRR10351810-CDS_Mapping_summary.txt \
$Sequence_length

java -cp ${script}/ Line_by_line_Total_dict_Get_exon_reads \
${file_output}/SRR10351811.sam \
SRR10351811 \
${file_output}/Tumor-SRR10351811-CDS_Mapping_summary.txt \
$Sequence_length

cat ${file_output}/Normal-SRR10351810-CDS_Mapping_summary.txt >> $summary_output/Total_WES_BWA_CDS_Glioma_normal.txt
cat ${file_output}/Tumor-SRR10351811-CDS_Mapping_summary.txt >>  $summary_output/Total_WES_BWA_CDS_Glioma_tumor.txt

##### Sequence Reads Mapping Quality #####

cd ${file_output}/

cat ${file_output}/SRR10351810.sam | cut -f5 > ${file_output}/SRR10351810-mapping_quality

Ngt30=$(cat SRR10351810-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Ngt60=$(cat SRR10351810-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ntotal=$(cat SRR10351810-mapping_quality| wc -l)
Nfra30=$(echo "$((Ngt30))/$((Ntotal))" | bc -l)
Nfra60=$(echo "$((Ngt60))/$((Ntotal))" | bc -l)

printf  "%s\t%4f\t%4f\n" "SRR10351810" "${Nfra30}" "${Nfra60}" >> $summary_output/WES_Total_Mapping_quality_Normal_Glioma.txt

cat ${file_output}/SRR10351811.sam | cut -f5 > ${file_output}/SRR10351811-mapping_quality

Tgt30=$(cat ${file_output}/SRR10351811-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Tgt60=$(cat ${file_output}/SRR10351811-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ttotal=$(cat ${file_output}/SRR10351811-mapping_quality| wc -l)
Tfra30=$(echo "$((Tgt30))/$((Ttotal))" | bc -l)
Tfra60=$(echo "$((Tgt60))/$((Ttotal))" | bc -l)

printf  "%s\t%4f\t%4f\n" "SRR10351811" "${Tfra30}" "${Tfra60}" >> $summary_output/WES_Total_Mapping_quality_Tumor_Glioma.txt

###### QC of Depth of Coverage ######

cd ${file_output}/

cat ${file_output}/SRR10351810_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${file_output}/SRR10351810_DepthofCoverage_Distribution.txt
cat ${file_output}/SRR10351811_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${file_output}/SRR10351811_DepthofCoverage_Distribution.txt

#cat SRR10351811_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
#cat SRR10351810_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
python ${script}/Randomness.py ${file_output}/SRR10351810_DepthofCoverage_Distribution.txt SRR10351810 Glioma Normal
python ${script}/Randomness.py ${file_output}/SRR10351811_DepthofCoverage_Distribution.txt SRR10351811 Glioma Tumor

cat ${file_output}/SRR10351810_randomness_summary.txt >> $summary_output/Total_WES_Glioma_normal_Randomness_Summary.txt
cat ${file_output}/SRR10351811_randomness_summary.txt >> $summary_output/Total_WES_Glioma_tumor_Randomness_Summary.txt

###### Callable bases ########


reference='/work/szlab/Lab_shared_PanCancer/source'
#Yuan_script='/work/szlab/Lab_shared_PanCancer/script'

cd ${file_output}

callable=$(cat i_03A6_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1"| wc -l)

printf  "%s\t%d\n" "i_03A6" "${callable}"  >> $summary_output/Total_WES_Callable_Glioma.txt


rm ${file_output}/SRR10351810.sam
rm ${file_output}/SRR10351811.sam
rm ${file_output}/SRR10351810_BWA_summarize
rm ${file_output}/SRR10351811_BWA_summarize
rm ${file_output}/SRR10351811_randomness_summary.txt
rm ${file_output}/SRR10351810_randomness_summary.txt
rm ${file_output}/Normal-SRR10351810-CDS_Mapping_summary.txt
rm ${file_output}/Tumor-SRR10351811-CDS_Mapping_summary.txt
rm ${file_output}/SRR10351810-mapping_quality ${file_output}/SRR10351811-mapping_quality
rm ${file_output}/SRR10351810_DepthofCoverage_Distribution.txt
${file_output}/SRR10351811_DepthofCoverage_Distribution.txt