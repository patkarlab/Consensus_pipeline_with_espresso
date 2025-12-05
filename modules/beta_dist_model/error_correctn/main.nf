#!/usr/bin/env nextflow

process ERROR_CORRECTN {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), path(annovar_snps)
		val(variant)
		path(matrix)
	output:
		tuple val (Sample), path("${Sample}_${variant}_error_corrected.csv")
	script:
	"""
	cal_beta_P_vals_vcf.py --matrix_file ${matrix} --input_csv_file ${annovar_snps} --output_csv_file ${Sample}_${variant}_error_corrected.csv
	"""
}