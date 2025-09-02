#!/usr/bin/env nextflow

process FILT3R {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*filt3r_json.csv'
	input:
		tuple val (Sample), file (trimmed_read1), file (trimmed_read2)
		filt3r_reference
	output:
		tuple val (Sample), file("*_filt3r.vcf"), file("*_filt3r_json.csv"), file("*_filt3r_out.csv")
	script:
	"""
	filt3r -k 12 --ref ${filt3r_reference} --sequences ${trimmed_read1},${trimmed_read2} --nb-threads ${task.cpus} --vcf --out ${Sample}_filt3r.json
	python3 convert_json_to_csv.py ${Sample}_filt3r.json ${Sample}_filt3r_json.csv
	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}_filt3r.vcf --outfile ${Sample}.filt3r.avinput --withzyg --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.filt3r.avinput --out ${Sample}.filt3r --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	# Check if the multianno file is empty
	if [[ ! -s ${Sample}.filt3r.hg19_multianno.csv ]]; then
		touch ${Sample}.filt3r__final.csv
		touch ${Sample}_filt3r_json_filtered.csv
		touch ${Sample}_filt3r_out.csv
	else
		somaticseqoutput-format_filt3r.py ${Sample}.filt3r.hg19_multianno.csv ${Sample}.filt3r__final.csv
		filter_json.py ${Sample}_filt3r_json.csv ${Sample}_filt3r_json_filtered.csv
		merge_filt3r_csvs.py ${Sample}.filt3r__final.csv ${Sample}_filt3r_json_filtered.csv ${Sample}_filt3r_out.csv
	fi
	"""
}

process MINIMAP_GETITD {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val (Sample), file (read1), file (read2)
		minimap_getitd_reference
	output:
		path "${Sample}_getitd"
	script:
	"""
	minimap2 -ax sr ${pminimap_getitd_reference} ${read1} ${read2} > ${Sample}.sam
	${params.samtools} view -b -h ${Sample}.sam -o ${Sample}.bam
	${params.samtools} sort ${Sample}.bam -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}

process GETITD {
	tag "${Sample}"
	publishDir "${params.outdir}/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val (Sample), file(finalBams), file(finalBamBai)
	output:
		path "${Sample}_getitd"
	script:
	"""
	${params.samtools} view ${finalBams} -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}

process COMBINE_FILT3R_GETITD {
	publishDir "{params.outdir}/${Sample}/", mode: 'copy'	
	input:
	tuple val (Sample), path (getitd_dir), file (combined_excel)  

	output:
	tuple val(Sample), file("${Sample}_Final_merged.xlsx")

	script:
	"""
	if [ -f "${getitd_dir}/itds_collapsed-is-same_is-similar_is-close_is-same_trailing.tsv" ]; then
		cp "${getitd_dir}/itds_collapsed-is-same_is-similar_is-close_is-same_trailing.tsv" "${Sample}_getitd.tsv"
	else
		getitd_gen.py -o "${Sample}_getitd.tsv"
	fi

	combine_all.py --getitd_out ${Sample}_getitd.tsv --merged_out ${combined_excel} --output ${Sample}_Final_merged.xlsx
	"""
}