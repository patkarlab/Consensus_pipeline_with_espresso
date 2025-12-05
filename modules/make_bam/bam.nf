#!/usr/bin/nextflow
// file paths
genome_file = file("${params.genome}", checkIfExists: true)
adaptors_file = file("${params.adaptors}", checkIfExists:true)

process TRIM {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val(Sample), file(read1), file(read2)
		path(Adapt)
	output:
		tuple val(Sample), file("${Sample}_trim_R1.fastq"), file("${Sample}_trim_R2.fastq")
	script:
	"""
	${params.fastp} -i ${read1} -I ${read2} -o ${Sample}_trim_R1.fastq -O ${Sample}_trim_R2.fastq --adapter_fasta ${Adapt} -w ${task.cpus} -q 20
	"""
}

process MAPBAM {
	tag "${Sample}"
	maxForks 5
	label 'process_low'
	input:
		tuple val(Sample), file(trim1), file(trim2)
	output:
		tuple val(Sample), file ("${Sample}.bam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" \
	-M -t ${task.cpus} ${params.genome} ${trim1} ${trim2} | ${params.samtools} sort -@ ${task.cpus} -o ${Sample}.bam -
	"""
}

process INDEX {
	tag "${Sample}"
	maxForks 5
	label 'process_low'
	input:
		tuple val(Sample), file(BamFile)
	output:
		tuple val(Sample), file ("${Sample}.bam"), file ("${Sample}.bam.bai")
	script:
	"""
	samtools index ${BamFile} > ${Sample}.bam.bai
	"""
}

process MARK_DUPS {
	tag "${Sample}"
	maxForks 5
	label 'process_medium'
	conda '/home/miniconda3/envs/gatk_4'
	input:
		tuple val(Sample), file(sortd_bam)
		path (fasta)
	output:
		tuple val(Sample), file("${Sample}_markdups.bam"), file("${Sample}_marked_dup_metrics.txt")
	script:
	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}g -XX:-UsePerfData" \
	MarkDuplicatesSpark \
	-I ${sortd_bam} \
	-O ${Sample}_markdups.bam \
	-M ${Sample}_marked_dup_metrics.txt \
	--conf "spark.executor.cores=${task.cpus}" \
	--conf "spark.executor.memory=${task.memory.toGiga()}g" \
	--conf "spark.driver.memory=${task.memory.toGiga()}g" \
	--conf "spark.network.timeout=10000000" \
	--conf "spark.executor.heartbeatInterval=10000000" \
	--tmp-dir .
	"""
}

process BQSR {
	tag "${Sample}"
	input:
		tuple val(Sample), file(markdups_bam), file(markdups_metrics)
	output:
		tuple val(Sample), file("${Sample}_recal.table")
	script:
	"""
	${params.gatk} BaseRecalibrator \
		-I ${markdups_bam} \
		-R ${params.genome} \
		--known-sites ${params.site1} \
		--known-sites ${params.site2} \
		--bqsr-baq-gap-open-penalty 30.0 \
		-O ${Sample}_recal.table
	"""
}

process APPLY_BQSR {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_uncollapsed.bam'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_uncollapsed.bam.bai'
	tag "${Sample}"
	input:
		tuple val(Sample), file(markdups_bam), file(markdups_metrics), file(recal_table)
	output:
		tuple val(Sample), file("${Sample}_uncollapsed.bam"), file("${Sample}_uncollapsed.bam.bai")
	script:
	"""
	${params.gatk} ApplyBQSR \
		-R ${params.genome} \
		-I ${markdups_bam} \
		--bqsr-recal-file ${recal_table} \
		-O ${Sample}_uncollapsed.bam

	mv ${Sample}_uncollapsed.bai ${Sample}_uncollapsed.bam.bai
	"""
}

process ALIGNMENT_METRICS {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_alignment_summary_metrics.txt'
	tag "${Sample}"
	input:
		tuple val(Sample), file(final_bam), file(final_bam_bai)
	output:
		tuple val(Sample), file("${Sample}_alignment_summary_metrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectAlignmentSummaryMetrics \
		R=${params.genome} \
		I=${final_bam} \
		O=${Sample}_alignment_summary_metrics.txt
	"""
}

process INSERT_SIZE_METRICS {
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_insert_size_metrics.txt'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_insert_size_metrics.pdf'
	tag "${Sample}"
	input:
		tuple val(Sample), file(final_bam), file(final_bam_bai)
	output:
		tuple val(Sample), file("${Sample}_insert_size_metrics.txt"), file("${Sample}_insert_size_metrics.pdf")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectInsertSizeMetrics \
		I=${final_bam} \
		O=${Sample}_insert_size_metrics.txt \
		H=${Sample}_insert_size_metrics.pdf \
		HISTOGRAM_WIDTH=500 \
		TMP_DIR=${Sample}_tmp
	"""
}

workflow FASTQTOBAM {
	take:
		samples_ch
	main:
		TRIM(samples_ch, adaptors_file)
		MAPBAM(TRIM.out)
		INDEX(MAPBAM.out)
		// MARK_DUPS(MAPBAM.out, genome_file)
		// BQSR(MARK_DUPS.out)
		// APPLY_BQSR(MARK_DUPS.out.join(BQSR.out))
		ALIGNMENT_METRICS(INDEX.out)
		INSERT_SIZE_METRICS(INDEX.out)
	emit:
		final_bams_ch = INDEX.out
		trim_ch = TRIM.out
}