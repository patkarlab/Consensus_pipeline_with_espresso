#!/usr/bin/nextflow

process somaticSeq_run {
	tag "${Sample}"
	label 'process_low'
	conda '/home/miniconda3/envs/new_base/'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.vcf'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_ErrorCorrectd.xlsx'
	input:
		tuple val (Sample), file (Coverage), file (mutectVcf), file (vardictVcf), file (varscanVcf), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("${Sample}.ErrorCorrectd.vcf"), emit: correctd_vcf_files
		tuple val (Sample), file ("${Sample}_ErrorCorrectd.xlsx"), emit: correctd_excels
	script:
	"""
	${params.vcf_filter} ${mutectVcf} ${Sample}.mutect2_sort.vcf
	${params.vcf_filter} ${vardictVcf} ${Sample}.vardict_sort.vcf
	${params.vcf_filter} ${varscanVcf} ${Sample}.varscan_sort.vcf
	echo "inside the somaticSeq"
	
	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads ${task.cpus} --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --mutect2-vcf ${Sample}.mutect2_sort.vcf --vardict-vcf ${Sample}.vardict_sort.vcf --varscan-vcf ${Sample}.varscan_sort.vcf --sample-name ${Sample}

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
	final_output.py ${Sample}_ErrorCorrectd.xlsx ${Sample}.ErrorCorrectd.csv ${Coverage}
	"""
}

process somaticSeq_run_uncoll {
	tag "${Sample}"
	label 'process_low'
	conda '/home/miniconda3/envs/new_base/'
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_uncollapsed.xlsx'
	input:
		tuple val (Sample), file (Coverage_uncollaps), file (mutectVcf), file (vardictVcf), file (varscanVcf), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.vcf_filter} ${mutectVcf} ${Sample}.mutect2_sort.vcf
	${params.vcf_filter} ${vardictVcf} ${Sample}.vardict_sort.vcf
	${params.vcf_filter} ${varscanVcf} ${Sample}.varscan_sort.vcf

	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads ${task.cpus} --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --mutect2-vcf ${Sample}.mutect2_sort.vcf --vardict-vcf ${Sample}.vardict_sort.vcf --varscan-vcf ${Sample}.varscan_sort.vcf --sample-name ${Sample}

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
	final_output.py ${Sample}_uncollapsed.xlsx ${Sample}.ErrorCorrectd.csv ${Coverage_uncollaps}
	"""
}