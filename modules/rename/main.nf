#!/usr/bin/nextflow

process RENAME_SMMIPS {
	tag "${Sample}"
	label 'process_medium'
	maxForks 5
	input:
		tuple val(Sample), path(assembled_fastq)
	output:
		tuple val(Sample), path("${Sample}_renamed.fastq.gz")
	script:
	"""
	CURRDIR=\$(pwd)
	TMPDIR=\$(mktemp -d)
	cp ${assembled_fastq} \$TMPDIR
	cd \$TMPDIR
	#gzip -dc ${assembled_fastq} > ${Sample}.fastq
	rename_smips_reads.pl ${assembled_fastq} ${Sample}_renamed.fastq
	#rm ${Sample}.fastq
	gzip ${Sample}_renamed.fastq
	
	mv ${Sample}_renamed.fastq.gz \$CURRDIR/
	"""
}
