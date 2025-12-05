#!/usr/bin/nextflow

process ALIGNMENT {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(renamed_fastq)
		path (GenFile)
		path (GenDir)
	output:
		tuple val(Sample), path("${Sample}_sortd.bam"), path("${Sample}_sortd.bam.bai"), emit: bam_out_ch
		// tuple val(Sample), path("${Sample}_chr13_sorted.bam"), emit: chr13_out_ch
	script:
	"""
	bwa mem -t ${task.cpus} ${GenFile} ${renamed_fastq} | samtools sort -@ ${task.cpus} -o ${Sample}_sortd.bam -
	samtools index ${Sample}_sortd.bam > ${Sample}_sortd.bam.bai

	#samtools view -@ ${task.cpus} ${Sample}_sortd.bam -b -h chr13 > ${Sample}_chr13.bam
	#samtools sort -@ ${task.cpus} ${Sample}_chr13.bam -o ${Sample}_chr13_sorted.bam
	"""
}