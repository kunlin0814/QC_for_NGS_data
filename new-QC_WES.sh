#!/bin/bash
#PBS -q batch
#PBS -N QC-WES-i_03A6
#PBS -l nodes=1:ppn=1
#PBS -l walltime=50:00:00
#PBS -l mem=55gb


# WGS 100 bps 
# WES 76 bps

result='/scratch/kh31516/Pan_cancer/glioma/results/WES/i_03A6' #
script='/work/szlab/kh31516_Lab_Share_script'
output_folder='/scratch/kh31516/Pan_cancer'
Sequence_length=76

module load SAMtools/1.9-foss-2016b
module load Anaconda3/2018.12
source activate py35

#
###### Analyze the BWA mapping result ######
###### BWA Mapping QC Columns ######
## ID, Total_pairs, Uniquely_mapped_rate, Repeatedly_mapped_rate, 1read_mapped_rate, Incorrectly_mapped_rate, Unmapped_rate, 
## Uniquely_mapped_count, Repeatedly_mapped_count, 1read_mapped_count, Incorrectly_mapped_count, Unmapped_count		
######

cd ${result}/

samtools view SRR10351810.bam > SRR10351810.sam # normal
samtools view SRR10351811.bam > SRR10351811.sam # tumor

java -cp ${script}/ Summarized_BWA \
${result}/SRR10351810.sam \
SRR10351810 \
SRR10351810_BWA_summarize

java -cp ${script}/ Summarized_BWA \
${result}/SRR10351811.sam \
SRR10351811 \
SRR10351811_BWA_summarize

cat ${result}/SRR10351810_BWA_summarize >> $output_folder/Total_WES_Normal_Glioma_BWA_summarization.txt
cat ${result}/SRR10351811_BWA_summarize >> $output_folder/Total_WES_Tumor_Glioma_BWA_summarization.txt


##### Sequence Reads Mapping Quality #####

cd ${result}/

cat ${result}/SRR10351810.sam | cut -f5 > ${result}/SRR10351810-mapping_quality

Ngt30=$(cat SRR10351810-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Ngt60=$(cat SRR10351810-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ntotal=$(cat SRR10351810-mapping_quality| wc -l)
Nfra30=$(echo "$((Ngt30))/$((Ntotal))" | bc -l)
Nfra60=$(echo "$((Ngt60))/$((Ntotal))" | bc -l)

printf  "%s\t%4f\t%4f\n" "SRR10351810" "${Nfra30}" "${Nfra60}" >> $output_folder/WES_Total_Mapping_quality_Normal_Glioma.txt

cat ${result}/SRR10351811.sam | cut -f5 > ${result}/SRR10351811-mapping_quality

Tgt30=$(cat ${result}/SRR10351811-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Tgt60=$(cat ${result}/SRR10351811-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ttotal=$(cat ${result}/SRR10351811-mapping_quality| wc -l)
Tfra30=$(echo "$((Tgt30))/$((Ttotal))" | bc -l)
Tfra60=$(echo "$((Tgt60))/$((Ttotal))" | bc -l)

printf  "%s\t%4f\t%4f\n" "SRR10351811" "${Tfra30}" "${Tfra60}" >> $output_folder/WES_Total_Mapping_quality_Tumor_Glioma.txt

###### CDS Regions Mapping Rates ######
#### The column of CDS mapping #####
#### file_name, Total_reads, Total_uniq, uniq_mapped_rate, Total_read_pairs, Uniq_Exonic_region_count, Uniq_Exonic_region_paris_rates

cd ${result}/

java -cp ${script}/ Total_dict_Get_exon_reads \
${result}/SRR10351810.sam \
SRR10351810 \
${result}/Normal-SRR10351810-CDS_Mapping_summary.txt \
${Sequence_length}

java -cp ${script}/ Total_dict_Get_exon_reads \
${result}/SRR10351811.sam \
SRR10351811 \
${result}/Tumor-SRR10351811-CDS_Mapping_summary.txt \
${Sequence_length}

cat ${result}/Normal-SRR10351810-CDS_Mapping_summary.txt >> $output_folder/Total_WES_CDS_Glioma_normal.txt
cat ${result}/Tumor-SRR10351811-CDS_Mapping_summary.txt >>  $output_folder/Total_WES_CDS_Glioma_tumor.txt


###### QC of Depth of Coverage ######
###### The Column Depth of Coverage #######
## file_name, average, standard_devitation, rmse, sumOfSqerror, Cancer_type, Status
#####################################

cd ${result}/

cat ${result}/SRR10351810_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${result}/SRR10351810_DepthofCoverage_Distribution.txt
cat ${result}/SRR10351811_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${result}/SRR10351811_DepthofCoverage_Distribution.txt

#cat SRR10351811_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
#cat SRR10351810_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
python ${script}/Randomness.py ${result}/SRR10351810_DepthofCoverage_Distribution.txt SRR10351810 Glioma Normal
python ${script}/Randomness.py ${result}/SRR10351811_DepthofCoverage_Distribution.txt SRR10351811 Glioma Tumor

cat ${result}/SRR10351810_randomness_summary.txt >> $output_folder/Total_WES_Glioma_normal_Randomness_Summary.txt
cat ${result}/SRR10351811_randomness_summary.txt >> $output_folder/Total_WES_Glioma_tumor_Randomness_Summary.txt



###### Callable bases ########


reference='/work/szlab/Lab_shared_PanCancer/source'
Yuan_script='/work/szlab/Lab_shared_PanCancer/script'

cd ${result}

callable=$(cat i_03A6_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1"| wc -l)

printf  "%s\t%d\n" "i_03A6" "${callable}"  >> $output_folder/Total_WES_Callable_Glioma.txt


rm ${result}/SRR10351810.sam
rm ${result}/SRR10351811.sam
rm ${result}/SRR10351810_BWA_summarize
rm ${result}/SRR10351811_BWA_summarize
rm ${result}/SRR10351811_randomness_summary.txt
rm ${result}/SRR10351810_randomness_summary.txt
rm ${result}/Normal-SRR10351810-CDS_Mapping_summary.txt
rm ${result}/Tumor-SRR10351811-CDS_Mapping_summary.txt
rm ${result}/SRR10351810-mapping_quality
rm ${result}/SRR10351811-mapping_quality