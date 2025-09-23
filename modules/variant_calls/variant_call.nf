#!/usr/bin/nextflow

process MUTECT2 {
	tag "${Sample}"
	maxForks 10
	label 'process_low'
	conda '/home/miniconda3/envs/gatk_4'
	// publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_mutect2.vcf'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}_mutect2.vcf")
	script:
	"""
	#${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} Mutect2 -R ${params.genome} -I:tumor ${abra_bam} -O ${Sample}.mutect2.vcf -L ${params.bedfile}.bed -mbq 20
	#${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} FilterMutectCalls -V ${Sample}.mutect2.vcf --stats ${Sample}.mutect2.vcf.stats -O ${Sample}_filtered.vcf -R ${params.genome} --unique-alt-read-count 3 --min-median-base-quality 20 --min-median-mapping-quality 30

	#${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} Mutect2 -R ${params.genome} -I:tumor ${abra_bam} -O ${Sample}.mutect2.vcf -L ${params.bedfile}.bed -mbq 20 --panel-of-normals ${params.mutect2_pon}
	#${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} FilterMutectCalls -V ${Sample}.mutect2.vcf --stats ${Sample}.mutect2.vcf.stats -O ${Sample}_filtered.vcf -R ${params.genome} --unique-alt-read-count 3 --min-median-base-quality 20 --min-median-mapping-quality 30
	#${params.java_path}/java -Xmx80G -jar ${params.GATK42_path} SelectVariants -V ${Sample}_filtered.vcf --exclude-filtered -O ${Sample}_selected.vcf

	# gatk 4 
	gatk Mutect2 -R ${params.genome} -I ${abra_bam} -O ${Sample}_mutect2.vcf -L ${params.bedfile}.bed --max-reads-per-alignment-start 0 --af-of-alleles-not-in-resource 1e-6
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
	// publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.varscan.vcf'
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

process DEEPSOMATIC {
	tag "${Sample}"
	input:
		tuple val (Sample), file(finalBam), file (finalBamBai)
	output:
		tuple val(Sample), file("*${Sample}_DS.vcf")
	script:
	""" 
	mkdir output
	outpath=`realpath output`
	bam_path=`realpath ${finalBam} | awk 'BEGIN{OFS=FS="/"} {\$NF=""; print \$0}'`
	vcf_output=${Sample}_DS.vcf
	control_bam_path=`realpath ${params.control_bam}`
	echo \$bam_path \$outpath \$vcf_output ${finalBam} \${control_bam_path} ${params.genome} ${params.bedfile}.bed
	${params.deepsomatic} \$bam_path \$outpath \$vcf_output ${finalBam} \${control_bam_path} ${params.genome} ${params.bedfile}.bed
	"""
}

process ANNOVAR {
	tag "${Sample}"
	label 'process_low'
	conda '/home/miniconda3/envs/new_base/'
	// publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), file (Coverage), file (mutectVcf), file (vardictVcf), file (varscanVcf)
	output:
		tuple val (Sample), file ("${Sample}_ErrorCorrected.xlsx")
	script:
	"""
	annovar.sh ${mutectVcf} ${varscanVcf} ${vardictVcf}
	somaticseqoutput-format_v2_mutect2.py mutect2.somaticseq.hg19_multianno.csv mutect2_.csv 
	somaticseqoutput-format_v2_varscan.py varscan.somaticseq.hg19_multianno.csv varscan_.csv
	somaticseqoutput-format_v2_vardict.py vardict.somaticseq.hg19_multianno.csv vardict_.csv
	combine_callers.py ${Sample}_ErrorCorrected.xlsx mutect2_.csv varscan_.csv vardict_.csv ${Coverage}
	"""
}

process ANNOVAR_UNCOLLAPSE {
	tag "${Sample}"
	label 'process_low'
	conda '/home/miniconda3/envs/new_base/'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), file (Coverage), file (mutectVcf), file (vardictVcf), file (varscanVcf)
	output:
		tuple val (Sample), file ("${Sample}_uncollaps.xlsx")
	script:
	"""
	annovar.sh ${mutectVcf} ${varscanVcf} ${vardictVcf}
	somaticseqoutput-format_v2_mutect2.py mutect2.somaticseq.hg19_multianno.csv mutect2_.csv 
	somaticseqoutput-format_v2_varscan.py varscan.somaticseq.hg19_multianno.csv varscan_.csv
	somaticseqoutput-format_v2_vardict.py vardict.somaticseq.hg19_multianno.csv vardict_.csv
	combine_callers.py ${Sample}_uncollaps.xlsx mutect2_.csv varscan_.csv vardict_.csv ${Coverage}
	"""
}