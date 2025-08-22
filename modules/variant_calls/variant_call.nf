#!/usr/bin/nextflow

process MUTECT2 {
	tag "${Sample}"
	maxForks 10
	label 'process_low'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}_selected.vcf")
	script:
	"""
	#${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} Mutect2 -R ${params.genome} -I:tumor ${abra_bam} -O ${Sample}.mutect2.vcf -L ${params.bedfile}.bed -mbq 20
	#${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} FilterMutectCalls -V ${Sample}.mutect2.vcf --stats ${Sample}.mutect2.vcf.stats -O ${Sample}_filtered.vcf -R ${params.genome} --unique-alt-read-count 3 --min-median-base-quality 20 --min-median-mapping-quality 30

	${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} Mutect2 -R ${params.genome} -I:tumor ${abra_bam} -O ${Sample}.mutect2.vcf -L ${params.bedfile}.bed -mbq 20 --panel-of-normals ${params.mutect2_pon}
	${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} FilterMutectCalls -V ${Sample}.mutect2.vcf --stats ${Sample}.mutect2.vcf.stats -O ${Sample}_filtered.vcf -R ${params.genome} --unique-alt-read-count 3 --min-median-base-quality 20 --min-median-mapping-quality 30
	${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} SelectVariants -V ${Sample}_filtered.vcf --exclude-filtered -O ${Sample}_selected.vcf
	"""
}

process VARDICT {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -th ${task.cpus} -f 0.0001 -r 8 -N ${Sample} -b ${abra_bam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf
	"""
}

process VARSCAN {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -d 1000000 -A -a -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup

	if [ -s ${Sample}.mpileup ]; then
		${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	else 
		echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${Sample}.varscan_snp.vcf
	fi

	if [ -s ${Sample}.mpileup ]; then
		${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	else
		echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${Sample}.varscan_indel.vcf
	fi

	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}