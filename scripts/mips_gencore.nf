#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

process fastp {
	input:
		val (Sample)
	output:
		tuple val (Sample), file("*R1_umi.fq.gz"), file("**R2_umi.fq.gz")
	script:
	"""	
	${params.fastp} -i ${params.sequences}/${Sample}_*R1_*.fastq.gz -o ${Sample}_R1_umi.fq.gz -I ${params.sequences}/${Sample}_*R2_*.fastq.gz -O ${Sample}_R2_umi.fq.gz -U --umi_loc per_read --umi_len 4
	"""
}

process fastqc {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*'
	input:
		val Sample
	output:
		path "*"
	"""
	${params.fastqc} -o ./ -f fastq ${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz
	"""
}

process mapping_reads {
	input:
		tuple val (Sample), file(R1_umi_fastq), file(R2_umi_fastq)
	output:
		tuple val (Sample), file ("*.sam")
	script:
	"""
	bwa mem -t 32 -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -K 150000000 -Y ${params.genome} ${R1_umi_fastq} ${R2_umi_fastq} > ${Sample}.sam
	"""
}

process sam_conversion {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.uncollaps.bam*'
	input:
		tuple val (Sample), file (samfile)	
	output:
		tuple val(Sample), file ("*.uncollaps.bam"), file ("*.uncollaps.bam.bai")
	script:
	"""
	${params.samtools} view -bT ${params.genome} ${samfile} > ${Sample}.bam
	${params.samtools} sort ${Sample}.bam > ${Sample}.uncollaps.bam
	${params.samtools} index ${Sample}.uncollaps.bam > ${Sample}.uncollaps.bam.bai
	"""
}

process GencoreConsensus {	
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy'
	input:
		tuple val (Sample), file(uncollaps_bam), file(uncollaps_bam_bai)
	output:
		tuple val (Sample), file("*_consensus.bam"), file("*_consensus.bam.bai"), file("*.json"), file("*.html")
	script:
	"""
	${params.gencore} -i ${uncollaps_bam} -o ${Sample}_consensus.bam -r ${params.genome} -b ${params.bedfile}.bed -s 4
	${params.samtools} index ${Sample}_consensus.bam > ${Sample}_consensus.bam.bai
	"""
}

process ABRA2_realign {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.bam*'
	input:
		tuple val (Sample), file (consensus_bam), file (consensus_bai), file (json), file (html)
	output:
		tuple val (Sample), file ("*.ErrorCorrectd.bam"), file ("*.ErrorCorrectd.bam.bai")
	script:
	"""
	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${consensus_bam} --out ${Sample}.ErrorCorrectd.bam --ref ${params.genome} --threads 8 --targets ${params.bedfile}.bed --tmpdir ./ > abra.log
	${params.samtools} index ${Sample}.ErrorCorrectd.bam > ${Sample}.ErrorCorrectd.bam.bai
	"""
}

process coverage_mosdepth {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_umi_cov.regions.bed'
	input:
		tuple val (Sample), file (abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}_umi_cov.regions.bed")
	script:
	"""
	${params.mosdepth_script} ${abra_bam} ${Sample}_umi_cov ${params.bedfile}.bed
	#${params.mosdepth_script} ${abra_bam} ${Sample}_exon_umi_cov ${params.bedfile2}.bed
	"""
}

process mutect2_run {
	maxForks 10
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("mutect2_.csv")
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T MuTect2 -R ${params.genome} -I:tumor ${abra_bam} -o ${Sample}.mutect2.vcf --dbsnp ${params.site2} -L ${params.bedfile}.bed -nct 25 -mbq 20
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.mutect2.vcf --outfile ${Sample}.mutect2.avinput -allsample -withfreq --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.mutect2.avinput --out ${Sample}.mutect2.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	if [ -s ${Sample}.mutect2.somaticseq.hg19_multianno.csv ]; then
		python3 ${PWD}/scripts/somaticseqoutput-format_v2_mutect2.py ${Sample}.mutect2.somaticseq.hg19_multianno.csv mutect2_.csv
	else
		touch mutect2_.csv
	fi	
	"""
}

process vardict {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.0001 -r 8 -N ${Sample} -b ${abra_bam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf
	"""
}

process varscan {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -d 1000000 -A -a -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_snp.vcf

	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_indel.vcf

	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process somaticSeq_run {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.vcf'
	//publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.csv'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_ErrorCorrectd.xlsx'
	input:
		tuple val (Sample), file (Coverage), file (mutectVcf), file (vardictVcf), file (varscanVcf), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.vcf_filter} ${vardictVcf} ${Sample}.vardict_sort.vcf
	${params.vcf_filter} ${varscanVcf} ${Sample}.varscan_sort.vcf
	
	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --vardict-vcf ${Sample}.vardict_sort.vcf --varscan-vcf ${Sample}.varscan_sort.vcf --sample-name ${Sample}

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sSNV.vcf ${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_snv.vcf > ${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_snv.vcf.gz

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sINDEL.vcf ${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_indel.vcf > ${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${Sample}.somaticseq/somaticseq_snv.vcf.gz ${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.ErrorCorrectd.vcf

	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.ErrorCorrectd.vcf --outfile ${Sample}.ErrorCorrectd.avinput -allsample -withfreq --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.ErrorCorrectd.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	if [ -s ${Sample}.somaticseq.hg19_multianno.csv ]; then
		python3 ${params.format_somaticseq_script} ${Sample}.somaticseq.hg19_multianno.csv ${Sample}.ErrorCorrectd.csv
	else
		touch ${Sample}.ErrorCorrectd.csv
	fi
	python3 ${PWD}/scripts/final_output.py ${Sample}_ErrorCorrectd.xlsx ${Sample}.ErrorCorrectd.csv ${Coverage} ${mutectVcf}
	"""
}

workflow MIPS_MRD_Gencore {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.set { samples_ch }

	main:
		fastp(samples_ch)
		fastqc(samples_ch)
		mapping_reads(fastp.out)
		sam_conversion(mapping_reads.out)
		GencoreConsensus(sam_conversion.out)
		ABRA2_realign(GencoreConsensus.out)
		coverage_mosdepth(ABRA2_realign.out)
		mutect2_run(ABRA2_realign.out)
		vardict(ABRA2_realign.out)
		varscan(ABRA2_realign.out)
		somaticSeq_run(coverage_mosdepth.out.join(mutect2_run.out.join(vardict.out.join(varscan.out.join(ABRA2_realign.out)))))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}