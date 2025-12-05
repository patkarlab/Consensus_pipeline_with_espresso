#!/usr/bin/nextflow

process SAMTOOLS {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(cons_sam)
	output:
		tuple val(Sample), path("*.bam"), path("*.bam.bai")
	script:
	"""
	for sam_file in ${cons_sam}
	do
		file_name=\$( basename \${sam_file} .sam)
		samtools view -@ ${task.cpus} -u -b -S \${sam_file} | samtools sort -@ ${task.cpus} '-' -o \${file_name}.bam
		samtools index -@ ${task.cpus} \${file_name}.bam
	done
	"""
}