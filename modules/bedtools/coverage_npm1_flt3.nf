#!/usr/bin/env nextflow

process COVERAGE {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai)
		bedfile
	output:
		tuple val (Sample), file ("${Sample}.counts.bed")
	script:
	"""
	${params.bedtools} coverage -counts -a ${bedfile} -b ${finalBams} > ${Sample}.counts.bed
	"""
}
