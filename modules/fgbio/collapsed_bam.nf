#!/usr/bin/nextflow

process FastqToBam {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file(R1_trim_fastq), file(R2_trim_fastq)
	output:
		tuple val (Sample), file ("*.unmapped.bam")
	script:	
	"""
	java -Xmx${task.memory.toGiga()}g -jar "/home/programs/fgbio/fgbio-2.0.1.jar" --compression 1 --async-io FastqToBam \
	--input ${R1_trim_fastq} ${R2_trim_fastq} --read-structures 4M+T 4M+T \
	--sample ${Sample} --library ${Sample} --output ${Sample}.unmapped.bam
	"""
}

process MapBam {
	tag "${Sample}"
	label 'process_low'	
	input:
		tuple val (Sample), file(unmapped_bam)
	output:
		tuple val (Sample), file("*.mapped.bam")
	script:
	"""
	# Align the data with bwa and recover headers and tags from the unmapped BAM
	${params.samtools} fastq ${unmapped_bam} | bwa mem -t ${task.cpus} -p -K 150000000 -Y ${params.genome} - | \
	java -Xmx${task.memory.toGiga()}g -jar ${params.fgbio_path} --compression 1 --async-io ZipperBams \
	--unmapped ${unmapped_bam} --ref ${params.genome} --output ${Sample}.mapped.bam
	"""
}

process GroupReadsByUmi {
	tag "${Sample}"
	maxForks 5
	label 'process_medium'
	conda '/home/miniconda3/envs/new_base/'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: "${Sample}_family.png"
	input:
		tuple val (Sample), file(mapped_bam)
	output: 
		tuple val (Sample), file ("*.grouped.bam"), emit : grouped_bam_ch
		path ("${Sample}_family.png")
	script:
	"""
	java -Xmx${task.memory.toGiga()}g -jar ${params.fgbio_path} --compression 1 --async-io GroupReadsByUmi --input ${mapped_bam} --strategy Adjacency --edits 1 \
	--output ${Sample}.grouped.bam --family-size-histogram ${Sample}.tag-family-sizes.txt
	plot_family_sizes.py ${Sample}.tag-family-sizes.txt ${Sample}_family
	"""
}

process CallMolecularConsensusReads {
	tag "${Sample}"
	maxForks 5
	label 'process_medium'
	input:
		tuple val (Sample), file (grouped_bam)
	output:
		tuple val (Sample), file ("*.cons.unmapped.bam")
	script:
	"""
	# This step generates unmapped consensus reads from the grouped reads and immediately filters them
	java -Xmx${task.memory.toGiga()}g -jar ${params.fgbio_path} --compression 0 CallMolecularConsensusReads --input ${grouped_bam} \
	--output /dev/stdout --min-reads 4 --threads ${task.cpus} | java -Xmx${task.memory.toGiga()}g -jar ${params.fgbio_path} --compression 1 \
	FilterConsensusReads --input /dev/stdin --output ${Sample}.cons.unmapped.bam --ref ${params.genome} \
	--min-reads 4 --min-base-quality 20 --max-base-error-rate 0.25
	"""
}

process FilterConsBam {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file (cons_unmapped_bam)
	output:
		tuple val (Sample), file ("*.cons.filtered.bam"), file ("*.cons.filtered.bam.bai")
	script:
	"""
	${params.samtools} fastq ${cons_unmapped_bam} | bwa mem -t ${task.cpus} -p -K 150000000 -Y ${params.genome} - | \
	java -Xmx${task.memory.toGiga()}g -jar ${params.fgbio_path} --compression 0 --async-io ZipperBams --unmapped ${cons_unmapped_bam} \
	--ref ${params.genome} --tags-to-reverse Consensus --tags-to-revcomp Consensus | ${params.samtools} sort --threads ${task.cpus} -o ${Sample}.cons.filtered.bam
	${params.samtools} index ${Sample}.cons.filtered.bam > ${Sample}.cons.filtered.bam.bai
	"""
}

process SyntheticFastq {
	tag "${Sample}"
	input:
		tuple val (Sample), file (cons_filt_bam), file (cons_filt_bam_bai)
	output:
		tuple val (Sample), file ("*.synreads.bam"), file ("*.synreads.bam.bai")
	script:
	"""
	# Making a synthetic fastq and bam based on the number of reads in the sample bam file
	${params.PosControlScript} ${cons_filt_bam} ${Sample}_inreads
	${params.fastq_bam} ${Sample}_inreads		# This script will output a file named *_inreads.fxd_sorted.bam

	# Merging the Sample bam file with positive control bam file
	mv ${cons_filt_bam} ${Sample}.cons.filtered_merge.bam
	${params.samtools} merge -f ${Sample}.synreads.bam ${Sample}.cons.filtered_merge.bam ${Sample}_inreads.fxd_sorted.bam
	${params.samtools} index ${Sample}.synreads.bam > ${Sample}.synreads.bam.bai
	"""
}