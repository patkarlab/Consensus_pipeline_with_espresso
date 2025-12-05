#!/usr/bin/nextflow

process SAMTOFASTQ {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(bam_files), path(bai_files)
	output:
		tuple val(Sample), path("*.fastq")
	script:
	"""
	for bam_file in ${bam_files}
	do
		outfile_name=\$( echo \${bam_file} | sed 's:.sam.*::g' )	# Extracting the original sam file name
		gatk --java-options "-Xmx${task.memory.toGiga()}g" SamToFastq I=\${bam_file} F=\${outfile_name}.fastq
	done
	"""
}