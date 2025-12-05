#!/usr/bin/nextflow

process MPILEUP {
	tag "${Sample}"
	label 'process_medium'
	publishDir "${params.outdir}/${Sample}/bamfolder/", mode: 'copy', pattern: "*.bam*"
	input:
		tuple val(Sample), path(renamed_fastq)
		path (GenFile)
		path (GenDir)
	output:
		tuple val(Sample), path("*.mpileup"), emit: mpileup_ch
		tuple val(Sample), path("*.bam"), path("*.bai"), emit: bam_bai_ch
	script:
	"""
	for fastq_files in ${renamed_fastq}
	do
		file_name=\$( basename \${fastq_files} .fastq)
		bwa mem -t ${task.cpus} ${GenFile} \${fastq_files} | samtools sort -@ ${task.cpus} -o \${file_name}.bam -
		samtools index \${file_name}.bam
		samtools mpileup -F 0.2 -f ${GenFile} -d 100000 \${file_name}.bam > \${file_name}.mpileup
		#call_variants_from_mpileup.pl \${file_name}.mpileup '0' \${file_name}.vcf
	done
	"""
}
