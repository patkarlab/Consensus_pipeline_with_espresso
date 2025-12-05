#!/usr/bin/nextflow

process CALL_CONSENSUS {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(nosingles_sam)
	output:
		tuple val(Sample), path("*_con.sam")
	script:
	"""
	out_suffix="_con.sam"
	for sam_file in ${nosingles_sam}
	do
		file_name=\$( basename \${sam_file} .sam)
		output=\${file_name}\${out_suffix}

		fgbio -Xmx${task.memory.toGiga()}g --tmp-dir=. \
		SortBam -s TemplateCoordinate --input=\${sam_file} --output=\${file_name}_sortd.bam 

		fgbio -Xmx${task.memory.toGiga()}g --tmp-dir=. \
		CallMolecularConsensusReads --input=\${file_name}_sortd.bam --output=\${output} --error-rate-post-umi=30 --min-reads=2 --min-input-base-quality 20 \
		--threads ${task.cpus}
	done
	"""
}