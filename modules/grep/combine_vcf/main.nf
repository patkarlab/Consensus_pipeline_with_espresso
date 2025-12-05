#!/usr/bin/nextflow

process COMBINE_VCF {
	tag "${Sample}"
	label 'process_single'
	//publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.vcf.*'
	input:
		tuple val(Sample), path(forward_vcf), path(reverse_vcf)
	output:
		tuple val(Sample), path("${Sample}_combine_all_vcfs_condensed.vcf.point"), emit: vcf_point_ch 
		tuple val(Sample), path("${Sample}_combine_all_vcfs_condensed.vcf.indel"), emit: vcf_indel_ch
	script:
	"""
	grep -h chr ${forward_vcf} ${reverse_vcf} | sort -k 1,1 -k 2,2n -k 4,4 -k 5,5 > ${Sample}_combine_all_vcfs_sorted.vcf
	combine_final_vcfs.pl ${Sample}_combine_all_vcfs_sorted.vcf ${Sample}_combine_all_vcfs_condensed.vcf
	split_vcf_indel.pl ${Sample}_combine_all_vcfs_condensed.vcf
	"""
}
