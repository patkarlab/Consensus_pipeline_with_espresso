#!/usr/bin/nextflow

process VARSCAN {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
		path (genome_loc)
		path (bedfile)
	output:
		tuple val(Sample), file("${Sample}.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -d 1000000 -B -A -a -l ${bedfile} -f ${genome_loc} ${finalBam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-reads2 1 --min-avg-qual 5 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}
