#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process somaticSeq_run {
	tag "${Sample}"
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.vcf'
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
	
	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --mutect2-vcf ${Sample}.mutect2_sort.vcf --vardict-vcf ${Sample}.vardict_sort.vcf --varscan-vcf ${Sample}.varscan_sort.vcf --sample-name ${Sample}

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
	python3 ${PWD}/scripts/final_output.py ${Sample}_ErrorCorrectd.xlsx ${Sample}.ErrorCorrectd.csv ${Coverage}

	"""
}

process ERROR_CAL {
	// publishDir "$PWD/Final_Output/", mode: 'copy'
	input:
		file (ExcelFiles)
	// output:
	// 	file("hsmetrics_probewise.txt") 
	script:
	"""
	for i in ${ExcelFiles}
	do
		echo \${i} >> file_names.txt
	done
	# Extracting the file names for BNC
	bnc=\$( grep -i 'bnc' file_names.txt )
	# Extracting the file names for samples
	samples=\$( grep -vi 'bnc' file_names.txt )

	/home/pipelines/Consensus_pipeline_with_espresso/scripts/Error_correctn/error_correction.py --samples \${samples} --normals \${bnc}
	for i in `cat ${params.input}`
	do 
		if [ -s \${i}_*_mod.xlsx ]; then
			cp \${i}_*_mod.xlsx $PWD/Final_Output/\${i}/
		fi
	done
	"""
}