#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process FASTP {
	tag "${Sample}"
	conda '/home/arpit/miniconda3'
	input:
		tuple val (Sample), file(R1_trim_fastq), file(R2_trim_fastq)
	output:
		tuple val (Sample), file("${Sample}_R1_umi.fq.gz"), file("${Sample}_R2_umi.fq.gz")
	script:
	"""
	fastp -i ${R1_trim_fastq} -o ${Sample}_R1_umi.fq.gz -I ${R2_trim_fastq} -O ${Sample}_R2_umi.fq.gz \
	-g -W 5 -q 30 -u 40 -x -3 -l 75 -c -e 30 \
	--cut_front_mean_quality 1
	"""
}

process FASTP_UMI {
	tag "${Sample}"
	conda '/home/arpit/miniconda3'
	input:
		tuple val (Sample), file(R1_trim_fastq), file(R2_trim_fastq)
	output:
		tuple val (Sample), file("${Sample}_R1_umi.fq.gz"), file("${Sample}_R2_umi.fq.gz")
	script:
	"""
	fastp -i ${R1_trim_fastq} -o ${Sample}_R1_umi.fq.gz -I ${R2_trim_fastq} -O ${Sample}_R2_umi.fq.gz \
	-g -W 5 -q 30 -u 40 -x -3 -l 75 -c -e 30 \
	--cut_front_mean_quality 1 \
	-U --umi_loc per_read --umi_len 8 
	"""
}

