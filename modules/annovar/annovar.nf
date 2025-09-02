#!/usr/bin/env nextflow

process ANNOVAR {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*.varscan_.csv'
	input:
		tuple val (Sample), file(varscanVcf)
	output:
		tuple val (Sample), file ("*.hg19_multianno.csv"), file("*.varscan_.csv")
	script:
	"""
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${varscanVcf}  --outfile ${Sample}.varscan.avinput --withzyg --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.varscan.avinput --out ${Sample}.varscan --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	if [ -s ${Sample}.varscan.hg19_multianno.csv ]; then
		somaticseqoutput-format_v2_varscan.py ${Sample}.varscan.hg19_multianno.csv ${Sample}.varscan_.csv
	else
		touch ${Sample}.varscan_.csv
	fi
	sleep 5s
	"""
}
