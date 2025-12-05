#!/usr/bin/nextflow

process GETITD {
	tag "${Sample}"
	label 'process_low'
	//publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val(Sample), path(chr13_fastq)
		file (getitd_amplicon)
		file (getitd_amplicon_kayser)
	output:
		tuple val(Sample), path ("${Sample}_getitd")
	script:
	"""
	python ${params.get_itd_path}/getitd.py -reference ${getitd_amplicon} \
	-anno ${getitd_amplicon_kayser} -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
	-reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern ${task.cpus} ${Sample} ${Sample}_chr13.fastq
	"""
}