#!/usr/bin/nextflow

process COMBINE_CALLERS {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), file (varscan_multianno), file(varscan_csv), file(filt3r_vcf), file(filt3r_json), file(filt3r_out_csv), file(coverage)
	output:
		tuple val (Sample), file("${Sample}.xlsx")
	script:
	"""
	merge-tsv.py ${Sample}.xlsx ${varscan_csv} ${filt3r_out_csv} ${coverage}
	"""
}
