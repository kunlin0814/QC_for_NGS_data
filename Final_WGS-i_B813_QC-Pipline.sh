#!/bin/bash
#PBS -q batch
#PBS -N QC-WGS-i_B813
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
#PBS -l mem=60gb

# WGS 100 bps 
# WES 76 bps

### The script need sequence lengths for each bam or sam files ####

result='/scratch/kh31516/Pan_cancer/glioma/results/WGS/i_B813' #
script='/work/szlab/kh31516_Lab_Share_script'
output_folder='/scratch/kh31516/Pan_cancer/glioma/results'



module load SAMtools/1.9-foss-2016b
module load Anaconda3/2018.12
#source activate py35

# convert bam files into sam file


cd $result
#check1=$(ls | grep "SRR10351592.sam")
#check2=$(ls | grep "SRR10351593.sam")

## Accept the original bam file and convert into SAM file ###

samtools view $result/SRR10351592.bam > $result/SRR10351592.sam # normal
samtools view $result/SRR10351593.bam > $result/SRR10351593.sam # tumor


###### Analyze the BWA and CDS mapping result ######
cd ${result}/

java -cp ${script}/ Line_by_line_Total_dict_Get_exon_reads \
${result}/SRR10351592.sam \
SRR10351592 \
${result}/Normal-SRR10351592-CDS_Mapping_summary.txt \
100 # Sequence length

java -cp ${script}/ Line_by_line_Total_dict_Get_exon_reads \
${result}/SRR10351593.sam \
SRR10351593 \
${result}/Tumor-SRR10351593-CDS_Mapping_summary.txt \
100 # Sequence Length

cat ${result}/Normal-SRR10351592-CDS_Mapping_summary.txt >> $output_folder/WGS_Total_BWA_CDS_Glioma_WES_normal.txt
cat ${result}/Tumor-SRR10351593-CDS_Mapping_summary.txt >>  $output_folder/WGS_Total_BWA_CDS_Glioma_WES_tumor.txt


##### Sequence Reads Mapping Quality #####

cd ${result}/

cat ${result}/SRR10351592.sam | cut -f5 > ${result}/SRR10351592-mapping_quality

Ngt30=$(cat SRR10351592-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Ngt60=$(cat SRR10351592-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ntotal=$(cat SRR10351592-mapping_quality| wc -l)
Nfra30=$(echo "$((Ngt30))/$((Ntotal))" | bc -l)
Nfra60=$(echo "$((Ngt60))/$((Ntotal))" | bc -l)

printf  "%s\t%4f\t%4f\n" "SRR10351592" "${Nfra30}" "${Nfra60}" >> $output_folder/WGS_Total_Mapping_quality_Normal_Glioma

cat ${result}/SRR10351593.sam | cut -f5 > ${result}/SRR10351593-mapping_quality

Tgt30=$(cat ${result}/SRR10351593-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Tgt60=$(cat ${result}/SRR10351593-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ttotal=$(cat ${result}/SRR10351593-mapping_quality| wc -l)
Tfra30=$(echo "$((Tgt30))/$((Ttotal))" | bc -l)
Tfra60=$(echo "$((Tgt60))/$((Ttotal))" | bc -l)

printf  "%s\t%4f\t%4f\n" "SRR10351593" "${Tfra30}" "${Tfra60}" >> $output_folder/WGS_Total_Mapping_quality_Tumor_Glioma


###### QC of Depth of Coverage ######

cd ${result}/

cat ${result}/SRR10351592_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${result}/SRR10351592_DepthofCoverage_Distribution.txt
cat ${result}/SRR10351593_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${result}/SRR10351593_DepthofCoverage_Distribution.txt

cat SRR10351593_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
cat SRR10351592_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
python ${script}/Randomness.py ${result}/SRR10351592_DepthofCoverage_Distribution.txt SRR10351592 Glioma Normal
python ${script}/Randomness.py ${result}/SRR10351593_DepthofCoverage_Distribution.txt SRR10351593 Glioma Tumor

cat ${result}/SRR10351592_randomness_summary.txt >> $output_folder/WGS_Total_Glioma_normal_Randomness_Summary.txt
cat ${result}/SRR10351593_randomness_summary.txt >> $output_folder/WGS_Total_Glioma_tumor_Randomness_Summary.txt

###### Callable bases ######

reference='/work/szlab/Lab_shared_PanCancer/source'
Yuan_script='/work/szlab/Lab_shared_PanCancer/script'

cd ${result}


####### Annovar #######
# Extract PASS records from vcf
# awk '$7 == "PASS" {print $0}' ${result}/i_B813_rg_added_sorted_dedupped_removed.MuTect.vcf > ${result}/i_B813_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# 5 Steps filtering
#rep -w KEEP ${result}/i_B813_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,26,27,38,39 > ${result}/i_B813_PASS.stat
#python ${Yuan_script}/Filter_MutectStat_5steps.py ${result}/i_B813_PASS.stat ${result}/i_B813_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
#perl $reference/annovar_CanFam3.1.99.gtf/convert2annovar.pl -format vcf4old ${result}/i_B813_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut > ${result}/i_B813_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput

# annovar annotate
#perl $reference/annovar_CanFam3.1.99.gtf/annotate_variation.pl --buildver canFam3 ${result}/i_B813_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $reference/annovar_CanFam3.1.99.gtf

# add gene name
#python ${Yuan_script}/Add_GeneName_N_Signature.py ${result}/i_B813_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function $reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

callable=$(cat i_B813_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1"| wc -l)
 
printf  "%s\t%d\n" "i_B813" "${callable}"  >> $output_folder/WGS_Total_Callable_Glioma.txt

rm ${result}/SRR10351592.sam ${result}/SRR10351593.sam
rm ${result}/SRR10351593_randomness_summary.txt ${result}/SRR10351592_randomness_summary.txt
rm ${result}/Normal-SRR10351592-CDS_Mapping_summary.txt ${result}/Tumor-SRR10351593-CDS_Mapping_summary.txt
rm ${result}/SRR10351592_DepthofCoverage_Distribution.txt ${result}/SRR10351593_DepthofCoverage_Distribution.txt