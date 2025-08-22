#!/usr/bin/nextflow

process COVERAGE_BEDTOOLS {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_umi_cov.regions_bedtools.bed'
	input:
		tuple val (Sample), file (abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}_umi_cov.regions_bedtools.bed")
	script:
	"""
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${abra_bam} > "${Sample}_umi_cov.regions_bedtools.bed"
	"""	
}

process COVERAGE_BEDTOOLS_UNCOLL {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_uncollaps_bedtools.bed'
	input:
		tuple val (Sample), file (uncollaps_bam), file (uncollaps_bam_bai)
	output:
		tuple val (Sample), file ("${Sample}_uncollaps_bedtools.bed")
	script:
	"""
	${params.bedtools} coverage -counts -a ${params.bedfile}.bed -b ${uncollaps_bam} > "${Sample}_uncollaps_bedtools.bed"
	#${params.bedtools} coverage -counts -a ${params.bedfile2}.bed -b ${uncollaps_bam} > "${Sample}_exon_uncollaps_bedtools.bed"
	"""	
}