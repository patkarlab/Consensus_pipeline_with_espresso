#!/usr/bin/bash
# This script will take bam as input and give csv as output

sample_list="/new_disk/MIPS_UNIPATH_080825/Final_Output/npm1_pos.samples"
bam_location="/new_disk/MIPS_UNIPATH_080825/Final_Output"
bedfile="/home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal210125_sortd"
genome="/home/reference_genomes/hg19_broad/hg19_all.fasta"

for samples in `cat ${sample_list}`
do
	echo ${samples}
	${params.samtools} mpileup -F 0.2 -A -l ${bedfile}.bed -f ${genome} -d 100000 ${bam_location}/${samples}/${samples}.ErrorCorrectd.bam > ${bam_location}/${samples}/${Sample}.mpileup
done
