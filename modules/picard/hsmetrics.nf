#!/usr/bin/nextflow

process HSMETRICS {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_hsmetrics.txt'
	input:
		tuple val(Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${abra_bam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	"""	
}

process HSMETRICS_UNCOLL {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_uncoll_hsmetrics.txt'
	input:
		tuple val(Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*_uncoll_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${abra_bam} O= ${Sample}_uncoll_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	"""
}