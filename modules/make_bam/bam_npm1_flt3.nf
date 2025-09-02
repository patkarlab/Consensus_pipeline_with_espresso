#!/usr/bin/nextflow

process TRIM {
	tag "${Sample}"
	input:
		tuple val (Sample), file(read1), file(read2)
		illumina_adapters
		nextera_adapters
	output:
		tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
	script:
	"""	
	trimmomatic PE \
	${read1} ${read2} \
	-baseout ${Sample}.fq.gz \
	ILLUMINACLIP:${illumina_adapters}:2:30:10:2:keepBothReads \
	ILLUMINACLIP:${nextera_adapters}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40	
	sleep 5s
	"""
}

process MAP {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val (Sample), file (trimmed_read1), file (trimmed_read2)
		genome_loc
	output:
		tuple val(Sample), file ("*.bam"), file ("*.bam.bai")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:300" -t ${task.cpus} ${genome_loc} ${trimmed_read1} ${trimmed_read2} | ${params.samtools} sort -@ ${task.cpus} -o ${Sample}.bam -
	${params.samtools} index ${Sample}.bam > ${Sample}.bam.bai
	"""
}
