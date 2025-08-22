#!/usr/bin/nextflow

process GETITD {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		path "*_getitd"
	script:
	"""
	${params.samtools} sort ${abra_bam} -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}