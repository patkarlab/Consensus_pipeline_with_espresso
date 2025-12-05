#!/usr/bin/nextflow

process PEAR_SELF_ASSEMBLE {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), path(read1), path(read2)
	output:
		tuple val(Sample), path("${Sample}.assembled.fastq.gz")
	script:
	"""
	pear -f ${read1} -r ${read2} -o ${Sample} -j ${task.cpus}
	gzip ${Sample}*.fastq
	"""
}
