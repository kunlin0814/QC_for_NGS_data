#!/bin/bash

### Please edit the following file first
# expected arguments
# sample_name
# bam file folder (source folder)
# sam_file_ouput_folder
# depth of coverage folders
# tumor and normal sample_id
# sqeuence length
# Mutect folder
# Cancer type
# bioproject

source_folder='/scratch/jlw87619/PanK9cancer/WGS_Osteosarcoma_PRJNA525883/results' #
sam_file_output='/scratch/kh31516/Pan_cancer/OSA_WGS/QC'
depthOfCoverage='/scratch/jlw87619/PanK9cancer/WGS_Osteosarcoma_PRJNA525883/results'
Mutect='/scratch/jlw87619/PanK9cancer/WGS_Osteosarcoma_PRJNA525883/results'
script_output='/scratch/kh31516/Pan_cancer/OSA_WGS/QC_script'
summary_output='/scratch/kh31516/Pan_cancer/OSA_WGS/'
Bioproject='PRJNA525883'
Sequence_length=151
Cancer_Type='OSA'

tumorNormalTable='/scratch/kh31516/Pan_cancer/OSA_WGS/source/WGSpairs_TgenOSA.txt'
fullPipelineScript='/home/kh31516/Test_qsub_V/WGS_WLOS_QC.sh'

sample_names=($(cat $tumorNormalTable | grep -v ^# | cut -f1))
normal_sample_ids=($(cat $tumorNormalTable | grep -v ^# | cut -f2))
tumor_sample_ids=($(cat $tumorNormalTable | grep -v ^# | cut -f3))
number_of_samples=${#sample_names[@]}
 
## Use bash for loop 
for (( i=0; i<$number_of_samples; i++ ));
do
	sample_name=${sample_names[$i]}
	normal_sample=${normal_sample_ids[$i]}
	tumor_sample=${tumor_sample_ids[$i]}
	job_name=${sample_name}_FullQC

	sbatch -J ${job_name} -D $script_output -o ${job_name}-%j.out -e ${job_name}-%j.err \
	--export=sample_name=$sample_name,output_dir=$sam_file_output,Normal_Sample=$normal_sample,Tumor_Sample=$tumor_sample,
	depth_of_coverage=$depthOfCoverage,summary_output=$summary_output,Mutect=$Mutect,Bioproject=$Bioproject,Sequence_length=$Sequence_length,\
	Cancer_Type=$Cancer_Type,source_folder=$source_folder \
	$fullPipelineScript
	

done
