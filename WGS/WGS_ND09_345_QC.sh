#!/bin/bash
#SBATCH --job-name=ND09_345_Summarized_BWA_WGS         # Job name (ND09_345_Summarized_BWA_WGS)
#SBATCH --partition=batch           # Queue name (batch)
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run in a single task on a single node
#SBATCH --cpus-per-task=1           # Number of CPU cores per task (4)
#SBATCH --mem=20G                   # Job memory limit (10 GB)
#SBATCH --time=40:00:00              # Time limit hrs:min:sec or days-hours:minutes:seconds
#SBATCH --output=ND09_345_Summarized_BWA_WGS.%j.out    # Standard output log
#SBATCH --error=ND09_345_Summarized_BWA_WGS.%j.err     # Standard error log


Cancer_Type='OM'
sample_name="ND09_345"
source_folder='/scratch/kh31516/Original_Melanoma/WGS/WGS/results/ND09_345' #
sam_file_output='/scratch/kh31516/Original_Melanoma/WGS/WGS/store/QC/ND09_345'
script='/home/kh31516/kh31516_Lab_Share_script'
summary_output='/scratch/kh31516/Original_Melanoma/WGS/WGS/store'
DepthOfCoverage='/scratch/kh31516/Original_Melanoma/WGS/WGS/store/DepthOfCoverage/ND09_345'
Mutect='/scratch/kh31516/Original_Melanoma/WGS/WGS/store/Mutect/ND09_345'
Sequence_length=104
Normal_Sample='SAMN07376268'
Tumor_Sample='SAMN07376261'
Bioproject='PRJNA389294'

module load SAMtools/1.9-GCC-8.3.0
ml Anaconda3/2020.02
#source activate py35
ml Java/1.8.0_241
mkdir -p $sam_file_output
## Convert the bam file into sam file

cd ${source_folder}/

#samtools view $source_folder/${Normal_Sample}.bam > $sam_file_output/${Normal_Sample}.sam # normal
#samtools view $source_folder/${Tumor_Sample}.bam > $sam_file_output/${Tumor_Sample}.sam # tumor

###### Analyze the BWA mapping and CDS region mapping result ######

cd ${sam_file_output}/

#String SamFile = args[0];
#String ID = args[1];
#String file_output = args[2];
#int read_length = Integer.parseInt(args[3]);
#String Cancer_Type = args[4];
#String Status = args[5];
<<'EOF'
java -cp ${script}/ SummarizeBwaWgs \
${sam_file_output}/${Normal_Sample}.sam \
${Normal_Sample} \
${sam_file_output}/Normal-${Normal_Sample}-WGS_Mapping_summary.txt \
$Sequence_length \
${Cancer_Type} \
"Normal"

java -cp ${script}/ SummarizeBwaWgs \
${sam_file_output}/${Tumor_Sample}.sam \
${Tumor_Sample} \
${sam_file_output}/Tumor-${Tumor_Sample}-WGS_Mapping_summary.txt \
$Sequence_length \
${Cancer_Type} \
"Tumor"

cat ${sam_file_output}/Normal-${Normal_Sample}-WGS_Mapping_summary.txt >> $summary_output/Total_WGS_BWA_${Bioproject}_${Cancer_Type}.txt
cat ${sam_file_output}/Tumor-${Tumor_Sample}-WGS_Mapping_summary.txt  >> $summary_output/Total_WGS_BWA_${Bioproject}_${Cancer_Type}.txt


##### Sequence Reads Mapping Quality #####

cd ${sam_file_output}/

cat ${sam_file_output}/${Normal_Sample}.sam | cut -f5 > ${sam_file_output}/${Normal_Sample}-mapping_quality

Ngt30=$(cat ${sam_file_output}/${Normal_Sample}-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Ngt60=$(cat ${sam_file_output}/${Normal_Sample}-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ntotal=$(cat ${sam_file_output}/${Normal_Sample}-mapping_quality| wc -l)
Nfra30=$(echo "$((Ngt30))/$((Ntotal))" | bc -l)
Nfra60=$(echo "$((Ngt60))/$((Ntotal))" | bc -l)

printf  "%s\t%4f\t%4f\t%s\t%s\n" "${Normal_Sample}" "${Nfra30}" "${Nfra60}" "${Cancer_Type}" "Normal" >> \
$summary_output/WGS_Total_Mapping_quality_${Bioproject}_${Cancer_Type}.txt

cat ${sam_file_output}/${Tumor_Sample}.sam | cut -f5 > ${sam_file_output}/${Tumor_Sample}-mapping_quality

Tgt30=$(cat ${sam_file_output}/${Tumor_Sample}-mapping_quality | awk '{if ($1>=30) {print $1}}' | wc -l)
Tgt60=$(cat ${sam_file_output}/${Tumor_Sample}-mapping_quality | awk '{if ($1>=60) {print $1}}' | wc -l)
Ttotal=$(cat ${sam_file_output}/${Tumor_Sample}-mapping_quality| wc -l)
Tfra30=$(echo "$((Tgt30))/$((Ttotal))" | bc -l)
Tfra60=$(echo "$((Tgt60))/$((Ttotal))" | bc -l)

printf   "%s\t%4f\t%4f\t%s\t%s\n" "${Tumor_Sample}" "${Tfra30}" "${Tfra60}" "${Cancer_Type}" "Tumor" >> \
$summary_output/WGS_Total_Mapping_quality_${Bioproject}_${Cancer_Type}.txt
###### QC of Depth of Coverage and Randomness ######

cd ${DepthOfCoverage}/

cat ${DepthOfCoverage}/${Normal_Sample}_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${sam_file_output}/${Normal_Sample}_DepthofCoverage_Distribution.txt
cat ${DepthOfCoverage}/${Tumor_Sample}_DepthofCoverage_CDS.bed | cut -f2 | sort -n | uniq -c | awk -F " " '{print $1,$2}' > ${sam_file_output}/${Tumor_Sample}_DepthofCoverage_Distribution.txt
EOF
cd ${DepthOfCoverage}/
#cat ${Tumor_Sample}_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc
#cat ${Normal_Sample}_DepthofCoverage_Distribution.txt | awk '{ print $1, $1 * $2 }'| awk '{print $2}' | paste -sd+ - | bc

# V2 python randomness scripts input parameters
# input_file = sys.argv[1]
# file_name = sys.argv[2]
# Cancer_type = sys.argv[3]
# Status = sys.argv[4]

cd ${sam_file_output}
# the header = ('file_name'+'\t'+'mean_coverage'+'\t'+'rmse'+'\t'+ 'Cancer_type'+'\t'+ 'Status'+'\n')
python ${script}/V2_randomness.py \
${sam_file_output}/${Normal_Sample}_DepthofCoverage_Distribution.txt \
${Normal_Sample} \
${Cancer_Type} \
Normal

python ${script}/V2_randomness.py \
${sam_file_output}/${Tumor_Sample}_DepthofCoverage_Distribution.txt \
${Tumor_Sample} \
${Cancer_Type} \
Tumor

tail -n +2 ${sam_file_output}/${Normal_Sample}_randomness_summary.txt >> $summary_output/Total_WGS_${Bioproject}_${Cancer_Type}_Randomness_Summary.txt
tail -n +2 ${sam_file_output}/${Tumor_Sample}_randomness_summary.txt >> $summary_output/Total_WGS_${Bioproject}_${Cancer_Type}_Randomness_Summary.txt

###### Callable bases ########


reference='/work/szlab/Lab_shared_PanCancer/source'
#Yuan_script='/work/szlab/Lab_shared_PanCancer/script'

cd ${sam_file_output}

callable=$(cat $Mutect/${sample_name}_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt | grep "^1"| wc -l)

printf  "%s\t%d\t%s\n" "${sample_name}" "${callable}" "${Cancer_Type}" >> $summary_output/Total_WGS_Callable_${Bioproject}_${Cancer_Type}.txt


: "
rm ${sam_file_output}/${Normal_Sample}.sam
rm ${sam_file_output}/${Tumor_Sample}.sam
rm ${sam_file_output}/Normal-${Normal_Sample}-CDS_Mapping_summary.txt
rm ${sam_file_output}/Tumor-${Tumor_Sample}-CDS_Mapping_summary.txt
rm ${sam_file_output}/${Normal_Sample}-mapping_quality ${sam_file_output}/${Tumor_Sample}-mapping_quality
rm ${DepthOfCoverage}/${Normal_Sample}_DepthofCoverage_Distribution.txt
rm ${DepthOfCoverage}/${Tumor_Sample}_DepthofCoverage_Distribution.txt
rm ${DepthOfCoverage}/${Tumor_Sample}_randomness_summary.txt
rm ${DepthOfCoverage}/${Normal_Sample}_randomness_summary.txt
"
