#!/usr/bin/nextflow

process ABRA2_realign {
	tag "${Sample}"
	label 'process_low'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.bam*'
	input:
		tuple val (Sample), file (filt_bam), file (filt_bai)
	output:
		tuple val (Sample), file ("*.ErrorCorrectd.bam"), file ("*.ErrorCorrectd.bam.bai")
	script:
	"""
	${params.java_path}/java -Xmx${task.memory.toGiga()}g -jar ${params.abra2_path}/abra2-2.23.jar --in ${filt_bam} --out ${Sample}.ErrorCorrectd.bam \
	--ref ${params.genome} --threads ${task.cpus} --targets ${params.bedfile}.bed --tmpdir ./ > abra.log
	${params.samtools} index ${Sample}.ErrorCorrectd.bam > ${Sample}.ErrorCorrectd.bam.bai
	"""
}