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

module load SAMtools/1.9-foss-2016b
module load Anaconda3/2018.12
source activate py35

#-Xmx2g

###### Analyze the BWA mapping result ######

cd ${result}/

samtools view SRR10351810.bam > SRR10351810.sam # normal
samtools view SRR10351811.bam > SRR10351811.sam # tumor

java -Xmx2g -cp ${script}/ Summarized_BWA \
${result}/SRR10351810.sam \
SRR10351810 \
SRR10351810_BWA_summarize

java -Xmx2g -cp ${script}/ Summarized_BWA \
${result}/SRR10351811.sam \
SRR10351811 \
SRR10351811_BWA_summarize

cat ${result}/SRR10351810_BWA_summarize >> /scratch/kh31516/Pan_cancer/glioma/results/Total_Normal_Glioma_BWA_summarization.txt
cat ${result}/SRR10351811_BWA_summarize >> /scratch/kh31516/Pan_cancer/glioma/results/Total_Tumor_Glioma_BWA_summarization.txt


##### Sequence Reads Mapping Quality #####

cd ${result}/

cat ${result}/SRR10351810.sam | cut -f5 > ${result}/SRR10351810-mapping_quality

Ngt30=$(cat SRR10351810-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Ngt60=$(cat SRR10351810-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ntotal=$(cat SRR10351810-mapping_quality| wc -l)
Nfra30=$(echo "$((Ngt30))/$((Ntotal))" | bc -l)
Nfra60=$(echo "$((Ngt60))/$((Ntotal))" | bc -l)

printf  "%s\t%4f\t%4f\n" "SRR10351810" "${Nfra30}" "${Nfra60}" >> /scratch/kh31516/Pan_cancer/glioma/results/Total_Mapping_quality_Normal_Glioma

cat ${result}/SRR10351811.sam | cut -f5 > ${result}/SRR10351811-mapping_quality

Tgt30=$(cat ${result}/SRR10351811-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Tgt60=$(cat ${result}/SRR10351811-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ttotal=$(cat ${result}/SRR10351811-mapping_quality| wc -l)
Tfra30=$(echo "$((Tgt30))/$((Ttotal))" | bc -l)
Tfra60=$(echo "$((Tgt60))/$((Ttotal))" | bc -l)

printf  "%s\t%4f\t%4f\n" "SRR10351811" "${Tfra30}" "${Tfra60}" >> /scratch/kh31516/Pan_cancer/glioma/results/Total_Mapping_quality_Tumor_Glioma

###### CDS Regions Mapping Rates ######
cd ${result}/

java -Xmx2g -cp ${script}/ Total_dict_Get_exon_reads \
${result}/SRR10351810.sam \
SRR10351810 \
${result}/Normal-SRR10351810-CDS_Mapping_summary.txt \
76

java -Xmx2g -cp ${script}/ Total_dict_Get_exon_reads \
${result}/SRR10351811.sam \
SRR10351811 \
${result}/Tumor-SRR10351811-CDS_Mapping_summary.txt \
76

cat ${result}/Normal-SRR10351810-CDS_Mapping_summary.txt >> /scratch/kh31516/Pan_cancer/glioma/results/Total_CDS_Glioma_WES_normal.txt
cat ${result}/Tumor-SRR10351811-CDS_Mapping_summary.txt >>  /scratch/kh31516/Pan_cancer/glioma/results/Total_CDS_Glioma_WES_tumor.txt

###### QC of Depth of Coverage ######

cd ${result}/

cat ${result}/SRR10351810_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${result}/SRR10351810_DepthofCoverage_Distribution.txt
cat ${result}/SRR10351811_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${result}/SRR10351811_DepthofCoverage_Distribution.txt

#cat SRR10351811_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
#cat SRR10351810_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
python ${script}/Randomness.py ${result}/SRR10351810_DepthofCoverage_Distribution.txt SRR10351810 Glioma Normal
python ${script}/Randomness.py ${result}/SRR10351811_DepthofCoverage_Distribution.txt SRR10351811 Glioma Tumor

cat ${result}/SRR10351810_randomness_summary.txt >> /scratch/kh31516/Pan_cancer/glioma/results/Total_Glioma_normal_Randomness_Summary.txt
cat ${result}/SRR10351811_randomness_summary.txt >> /scratch/kh31516/Pan_cancer/glioma/results/Total_Glioma_tumor_Randomness_Summary.txt

###### Callable bases ########


reference='/work/szlab/Lab_shared_PanCancer/source'
Yuan_script='/work/szlab/Lab_shared_PanCancer/script'

cd ${result}

####### Annovar #######
# Extract PASS records from vcf
awk '$7 == "PASS" {print $0}' ${result}/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf > ${result}/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# 5 Steps filtering
grep -w KEEP ${result}/i_03A6_rg_added_sorted_dedupped_removed.bam_call_stats.txt | cut -f1,2,26,27,38,39 > ${result}/i_03A6_PASS.stat
python ${Yuan_script}/Filter_MutectStat_5steps.py ${result}/i_03A6_PASS.stat ${result}/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS

# annovar input preparation
perl $reference/annovar_CanFam3.1.99.gtf/convert2annovar.pl -format vcf4old ${result}/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut > ${result}/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput

# annovar annotate
perl $reference/annovar_CanFam3.1.99.gtf/annotate_variation.pl --buildver canFam3 ${result}/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput $reference/annovar_CanFam3.1.99.gtf

# add gene name
python ${Yuan_script}/Add_GeneName_N_Signature.py ${result}/i_03A6_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function $reference/Canis_familiaris.CanFam3.1.99.chr.gtf_geneNamePair.txt

callable=$(cat i_03A6_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1"| wc -l)

printf  "%s\t%d\n" "i_03A6" "${callable}"  >> /scratch/kh31516/Pan_cancer/glioma/results/WES_Total_Callable_Glioma.txt


rm ${result}/SRR10351810.sam
rm ${result}/SRR10351811.sam
rm ${result}/SRR10351810_BWA_summarize
rm ${result}/SRR10351811_BWA_summarize
rm ${result}/SRR10351811_randomness_summary.txt
rm ${result}/SRR10351810_randomness_summary.txt
rm ${result}/Normal-SRR10351810-CDS_Mapping_summary.txt
rm ${result}/Tumor-SRR10351811-CDS_Mapping_summary.txt