#!/usr/bin/nextflow

process COMBINE_MIPS_COUNTS {
	tag "${Sample}"
	label 'process_single'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val(Sample), path(forward_mip_count), path(reverse_mip_count)
	output:
		tuple val(Sample), path("${Sample}_mip_counts.txt")
	script:
	"""
	cat ${forward_mip_count} ${reverse_mip_count} > ${Sample}_mip_counts.txt
	"""
}