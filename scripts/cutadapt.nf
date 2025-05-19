#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process CUTADAPT {
	tag "${Sample}"
	input:
		tuple val (Sample), file(trimmed_forward), file(trimmed_reverse)
	output:
		tuple val(Sample), file("${Sample}_trimmed_R1.fastq.gz"), file("${Sample}_trimmed_R2.fastq.gz")
	script:
	"""
	# The sequences to remove from the 5' end were obtained from supplementary data of Leukemia 33, 2241-2253 (2019)
	${params.cutadapt} \
	--minimum-length 40 -g GTAAAACGACGGCCAGT -G TAATACGACTCACTATAGGG \
	-o ${Sample}_trimmed_R1.fastq.gz -p ${Sample}_trimmed_R2.fastq.gz \
	${trimmed_forward} ${trimmed_reverse}
	"""
}

process CUTADAPT_FASTQ {
	tag "${Sample}"
	input:
		val (Sample)
	output:
		tuple val(Sample), file("${Sample}_trimmed_R1.fastq.gz"), file("${Sample}_trimmed_R2.fastq.gz")
	script:
	"""
	# The sequences to remove from the 5' end were obtained from supplementary data of Leukemia 33, 2241-2253 (2019)
	${params.cutadapt} \
	--minimum-length 40 -g GTAAAACGACGGCCAGT -G TAATACGACTCACTATAGGG \
	-o ${Sample}_trimmed_R1.fastq.gz -p ${Sample}_trimmed_R2.fastq.gz \
	${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz
	"""
}

process BAMTOFASTQ {
	tag "${Sample}"
	input:
		tuple val (Sample), file (BamFile), file (BamBaiFile)
	output:
		tuple val (Sample), file ("${Sample}_R1.fastq.gz"), file ("${Sample}_R2.fastq.gz")
	script:
	"""
	${params.samtools} fastq -@ 8 ${BamFile} -1 ${Sample}_R1.fastq.gz -2 ${Sample}_R2.fastq.gz -0 /dev/null -s /dev/null -n
	"""
}

process FastqToBam_ASSEMBLED {
	tag "${Sample}"
	input:
		tuple val (Sample), file (pairAssembled)
	output:
		tuple val (Sample), file ("*.unmapped.bam")
	script:	
	"""
	sed '1~4s/_//g' ${pairAssembled} > ${Sample}.fastq
	java -Xmx12g -jar ${params.fgbio_new} --compression 1 --async-io FastqToBam --input ${Sample}.fastq --extract-umis-from-read-names --sample ${Sample} --library ${Sample} --output ${Sample}.unmapped.bam
	"""
}