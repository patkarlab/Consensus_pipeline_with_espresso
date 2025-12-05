#!/usr/bin/nextflow

process SPLIT_MIPS {
	tag "${Sample}"
	label 'process_single'
	input:
		tuple val(Sample), path(forward_sam_list), path(reverse_sam_list)
		path(txtfile)
		val(sam_orientation)
		path(scripts)
	output:
		tuple val(Sample), path("*/*.nosingles.RX_withheader.sam"), emit: nosingles_sam_ch
		tuple val(Sample), path("${Sample}_${sam_orientation}_mip_counts.txt"), emit: mip_count_ch
	script:
	"""
	mips.sh ${Sample} ${txtfile} ${sam_orientation} ${Sample}_${sam_orientation}_mip_counts.txt
	"""
}