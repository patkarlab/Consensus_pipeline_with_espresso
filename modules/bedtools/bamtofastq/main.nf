#!/usr/bin/nextflow

process BAMTOFASTQ {
	tag "${Sample}"
	label 'process_single'
	input:
		tuple val(Sample), path(chr13_bam)
	output:
		tuple val(Sample), path("${Sample}_chr13.fastq")
	script:
	"""
	bedtools bamtofastq -i ${chr13_bam} -fq ${Sample}_chr13.fastq
	"""
}