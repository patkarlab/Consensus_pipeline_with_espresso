#!/usr/bin/nextflow

process SPLITSAM {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(bam_file), path(bai_file)
	output:
		tuple val(Sample), path ("${Sample}_chr*F.sam"), path ("${Sample}_chr*R.sam")
	script:
	"""
	for chr in \$(seq 1 22) X Y; do
		samtools view -@ ${task.cpus} -f 16 -F 2048 ${bam_file} "chr\${chr}" | \
		LC_ALL=C sort --parallel=${task.cpus} -S ${task.memory.toGiga()}G -k 4,4n > ${Sample}_"chr\${chr}"F.sam

		samtools view -@ ${task.cpus} -F 2064 ${bam_file} "chr\${chr}" | \
		LC_ALL=C sort --parallel=${task.cpus} -S ${task.memory.toGiga()}G -k 4,4n > ${Sample}_"chr\${chr}"R.sam
	done
	#split_sam.sh \${sam_file}F.sam \${sam_file}R.sam ${Sample}
	"""
}