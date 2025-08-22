#!/usr/bin/nextflow

process CNS_FILEGEN {
	tag "${Sample}"
	// publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.final.cns'
	input:
		tuple val (Sample), file (abra_bam), file (abra_bai)
	output:
		file ("*.final.cns")
	script:
	"""

	${params.samtools} mpileup -s -x -BQ 0 -q 1 -d 1000000 -A -a --skip-indels -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	# These parmeters were taken from the Espresso paper(Abelson et al., Sci. Adv. 2020; 6 : eabe3722)
	# A stringent (lower) p value removes all the snps hence, p-value is 1
	${params.java_path}/java -jar ${params.varscan_path} pileup2cns ${Sample}".mpileup" --variants SNP \
	--min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --strand-filter 0 > ${Sample}".cns"

	grep -v '^chrM' ${Sample}".cns" > ${Sample}".nochrM.cns"
	${params.filter_cns} ${Sample}".nochrM.cns" ${Sample}".final.cns"
	"""
}

process ESPRESSO {
	input:
		path (cns_file_list)
	script:
	"""
	sample_ids="${cns_file_list.join(' ')}"
	${params.umi_error_model} ${params.input} ${params.bedfile}.bed ${params.outdir}
	"""
}