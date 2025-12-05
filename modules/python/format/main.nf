#!/usr/bin/env nextflow

process FORMAT_ANNOVAR {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), val(vartype), file(anno_out)
		path (list_of_variants)
	output:
		tuple val (Sample), path("${Sample}_${vartype}.csv")
	script:
	"""
	format.py ${anno_out} ${Sample}.csv
	filter_variants.py --input_file ${Sample}.csv --variant_list ${list_of_variants} --output ${Sample}_${vartype}.csv
	"""
}

process FORMAT_ANNOVAR_SNP {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), val(vartype), file(anno_out)
	output:
		tuple val (Sample), path("${Sample}_${vartype}.csv")
	script:
	"""
	format.py ${anno_out} ${Sample}_${vartype}.csv
	"""
}