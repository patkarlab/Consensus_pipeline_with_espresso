#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

process ABRA2_realign {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.bam*'
	input:
		val Sample
	output:
		tuple val (Sample), file ("*.ErrorCorrectd.bam"), file ("*.ErrorCorrectd.bam.bai")
	script:
	"""
	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${params.sequences}/${Sample}.bam --out ${Sample}.ErrorCorrectd.bam --ref ${params.genome} --threads 8 --targets ${params.bedfile}.bed --tmpdir ./ > abra.log
	${params.samtools} index ${Sample}.ErrorCorrectd.bam > ${Sample}.ErrorCorrectd.bam.bai
	"""
}

process CNS_filegen {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.final.cns'
	input:
		tuple val (Sample), file (abra_bam), file (abra_bai)
	output:
		tuple val(Sample), file ("*.final.cns")
	script:
	"""
	#${params.samtools} mpileup -s -x -BQ 0 -q 1 -d 1000000 --skip-indels -f ${params.genome} ${abra_bam} > ${Sample}.mpileup

	#${params.samtools} mpileup -A -a --skip-indels -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	#${params.samtools} mpileup -d 1000000 -A -a --skip-indels -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup

	${params.samtools} mpileup -s -x -BQ 0 -q 1 -d 1000000 -A -a --skip-indels -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	# These parmeters were taken from the Espresso paper(Abelson et al., Sci. Adv. 2020; 6 : eabe3722)
	# A stringent (lower) p value removes all the snps hence, p-value is 1
	# ${params.java_path}/java -jar ${params.varscan_path} pileup2cns ${Sample}".mpileup" --variants SNP --min-coverage 10 --min-reads2 1 --min-avg-qual 30 --min-var-freq 0.0001 --p-value 1 --strand-filter 0 > ${Sample}".cns"

	#${params.java_path}/java -jar ${params.varscan_path} pileup2cns ${Sample}".mpileup" --variants SNP --min-coverage 2 --min-reads2 1 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1 --strand-filter 0 > ${Sample}".cns"

	${params.java_path}/java -jar ${params.varscan_path} pileup2cns ${Sample}".mpileup" --variants SNP --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --strand-filter 0 > ${Sample}".cns"

	grep -v '^chrM' ${Sample}".cns" > ${Sample}".nochrM.cns"
	${params.filter_cns} ${Sample}".nochrM.cns" ${Sample}".final.cns"
	"""
}

process espresso {
	input:
		val (tuple)
	script:	
	"""
	for i in `cat ${params.input}`
	do
		ln -s $PWD/Final_Output/\${i}/\${i}.final.cns ./
	done
	${params.umi_error_model} ${params.input} ${params.bedfile}.bed
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


process hsmetrics_run {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_hsmetrics.txt'
	input:
		tuple val(Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${abra_bam} O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
	"""	
}

process minimap_getitd {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_getitd'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		path "*_getitd"
	script:
	"""
	${params.samtools} sort ${abra_bam} -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam
	${params.bedtools} bamtofastq -i ${Sample}.chr13.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	"""
}

process mutect2_run {
	maxForks 10
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.mutect2.vcf")
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T MuTect2 -R ${params.genome} -I:tumor ${abra_bam} -o ${Sample}.mutect2.vcf --dbsnp ${params.site2} -L ${params.bedfile}.bed -nct 25 -mbq 20
	#${params.samtools} view -bs 40.1 ${abra_bam} > subsampled_01.bam
	#${params.samtools} index subsampled_01.bam
	#${params.mutect2} ${params.java_path} ${params.GATK38_path} ${params.genome} subsampled_01.bam ${Sample}.mutect2.vcf ${params.site2} ${params.bedfile}.bed
	"""
}

process vardict {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	#VarDict -G ${params.genome} -f 0.0001 -N ${Sample} -b ${abra_bam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf
	VarDict -G ${params.genome} -f 0.0001 -r 8 -N ${Sample} -b ${abra_bam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf
	"""
}

process lofreq {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val(Sample), file ("*.lofreq.filtered.vcf")
	script:
	"""
	${params.lofreq_path} viterbi -f ${params.genome} -o ${Sample}.lofreq.pre.bam ${abra_bam}
	${params.samtools} sort ${Sample}.lofreq.pre.bam > ${Sample}.lofreq.bam
	#${params.lofreq_path} call -b dynamic -C 2 -a 0.00005 -q 20 -Q 20 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam --no-default-filter
	${params.lofreq_path} call -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam --no-default-filter
	#${params.lofreq_path} filter --no-defaults -a 0.0001 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	${params.lofreq_path} filter -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}

process varscan {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.varscan.vcf")
	script:
	"""
	#${params.samtools} mpileup -d 1000000 -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	${params.samtools} mpileup -d 1000000 -A -a -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	#${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_snp.vcf

	#${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_indel.vcf

	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process strelka {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.strelka*.vcf'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.strelka*.vcf")
	script:
	"""
	${params.strelka_path}/configureStrelkaGermlineWorkflow.py --bam ${abra_bam} --referenceFasta ${params.genome} --callRegions ${params.bedfile}.bed.gz --targeted --runDir ./
	./runWorkflow.py -m local -j 20
	gunzip -f ./results/variants/variants.vcf.gz
	mv ./results/variants/variants.vcf ${Sample}.strelka.vcf

	${params.strelka_path}/configureStrelkaSomaticWorkflow.py --normalBam ${params.NA12878_bam} --tumorBam ${abra_bam} --referenceFasta ${params.genome} --callRegions ${params.bedfile}.bed.gz --targeted --runDir ./strelka-somatic
	./strelka-somatic/runWorkflow.py -m local -j 20
	${params.bcftools_path} concat -a strelka-somatic/results/variants/somatic.indels.vcf.gz strelka-somatic/results/variants/somatic.snvs.vcf.gz -o ./${Sample}.strelka-somatic.vcf
	"""
}

process mips_var_caller {
	//publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*mips_var_caller.csv'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		//tuple val (Sample), file ("*mips_var_caller.csv")
		tuple val (Sample), file ("*.mips.vcf")
	script:
	"""
	${params.samtools} mpileup -F 0.2 -A -l ${params.bedfile}.bed -f ${params.genome} -d 100000 ${abra_bam} > ${Sample}.mpileup
	#perl /home/pipelines/Consensus_pipeline_with_espresso/scripts/call_variants_from_mpileup.pl ${Sample}.mpileup 20 ${Sample}_mips.vcf
	/home/pipelines/Consensus_pipeline_with_espresso/scripts/mips_variant_call.sh ${Sample}.mpileup 20 ${Sample}.mips.vcf

	#perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}_mips.vcf --outfile temp.avinput --withzyg --includeinfo
	#perl ${params.annovarLatest_path}/table_annovar.pl temp.avinput --out ${Sample}.mips --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt
	#if [ ! -s ${Sample}.mips.hg19_multianno.csv ]
	#then
	#	touch ${Sample}_mips_var_caller.csv
	#else
	#	/home/pipelines/Consensus_pipeline_with_espresso/scripts/somaticseqoutput-format_v2_mips.py ${Sample}.mips.hg19_multianno.csv ${Sample}_mips_var_caller.csv
	#fi
	"""
}

process somaticSeq_run {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.vcf'
	//publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.csv'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_ErrorCorrectd.xlsx'
	input:
		tuple val (Sample), file (Coverage), file (mipsVcf), file (mutectVcf), file (vardictVcf), file (varscanVcf), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.vcf_filter} ${mutectVcf} ${Sample}.mutect2_sort.vcf
	${params.vcf_filter} ${vardictVcf} ${Sample}.vardict_sort.vcf
	${params.vcf_filter} ${varscanVcf} ${Sample}.varscan_sort.vcf
	${params.vcf_filter} ${mipsVcf} ${Sample}.mips_sort.vcf
	echo "inside the somaticSeq"
	
	python3 ${params.splitvcf_path} -infile ${Sample}.mips_sort.vcf -snv ${Sample}_mips_cnvs.vcf -indel ${Sample}_mips_indels.vcf
	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --mutect2-vcf ${Sample}.mutect2_sort.vcf --vardict-vcf ${Sample}.vardict_sort.vcf --varscan-vcf ${Sample}.varscan_sort.vcf --sample-name ${Sample} --arbitrary-snvs ${Sample}_mips_cnvs.vcf --arbitrary-indels ${Sample}_mips_indels.vcf

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sSNV.vcf ${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_snv.vcf > ${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_snv.vcf.gz

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sINDEL.vcf ${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_indel.vcf > ${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${Sample}.somaticseq/somaticseq_snv.vcf.gz ${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.ErrorCorrectd.vcf

	sed -i 's/##INFO=<ID=MVD0,Number=4,Type=Integer,Description="Calling decision of the 4 algorithms: MuTect, VarScan2, VarDict, SnvCaller_0">/##INFO=<ID=MVDS,Number=4,Type=String,Description="Calling decision of the 4 algorithms: MuTect, VarScan2, VarDict, MIPS_caller">/g' ${Sample}.ErrorCorrectd.vcf
	sed -i 's/MVD0/MVDS/g' ${Sample}.ErrorCorrectd.vcf

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

process somaticSeq_run_uncoll {
	//publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.UncollapsVariant.vcf'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_uncollapsed.xlsx'
	input:
		tuple val (Sample), file (Coverage_uncollaps), file (Coverage_uncollaps_exon), file (mutectVcf), file (vardictVcf), file (varscanVcf), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.vcf_filter} ${mutectVcf} ${Sample}.mutect2_sort.vcf
	${params.vcf_filter} ${vardictVcf} ${Sample}.vardict_sort.vcf
	${params.vcf_filter} ${varscanVcf} ${Sample}.varscan_sort.vcf

	#python3 ${params.splitvcf_path} -infile ${Sample}.mips_sort.vcf -snv ${Sample}_mips_cnvs.vcf -indel ${Sample}_mips_indels.vcf
	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --pass-threshold 0 --lowqual-threshold 0 --algorithm xgboost -minMQ 0 -minBQ 0 -mincaller 0 --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file ${finalBam} --mutect2-vcf ${Sample}.mutect2_sort.vcf --vardict-vcf ${Sample}.vardict_sort.vcf --varscan-vcf ${Sample}.varscan_sort.vcf --sample-name ${Sample}

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sSNV.vcf ${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_snv.vcf > ${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_snv.vcf.gz

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sINDEL.vcf ${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_indel.vcf > ${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${Sample}.somaticseq/somaticseq_snv.vcf.gz ${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.ErrorCorrectd.vcf

	#sed -i 's/##INFO=<ID=MVD0,Number=4,Type=Integer,Description="Calling decision of the 4 algorithms: MuTect, VarScan2, VarDict, SnvCaller_0">/##INFO=<ID=MVDS,Number=4,Type=String,Description="Calling decision of the 4 algorithms: MuTect, VarScan2, VarDict, MIPS_caller">/g' ${Sample}.ErrorCorrectd.vcf
	#sed -i 's/MVD0/MVDS/g' ${Sample}.ErrorCorrectd.vcf

	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.ErrorCorrectd.vcf --outfile ${Sample}.ErrorCorrectd.avinput -allsample -withfreq --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.ErrorCorrectd.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	if [ -s ${Sample}.somaticseq.hg19_multianno.csv ]; then
		python3 ${params.format_somaticseq_script} ${Sample}.somaticseq.hg19_multianno.csv ${Sample}.ErrorCorrectd.csv
	else
		touch ${Sample}.ErrorCorrectd.csv
	fi
	python3 ${PWD}/scripts/final_output.py ${Sample}_uncollapsed.xlsx ${Sample}.ErrorCorrectd.csv ${Coverage_uncollaps}
	"""
}

workflow NARASIMHA_MRD {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.set { samples_ch }

	main:
		ABRA2_realign(samples_ch)
		coverage_mosdepth(ABRA2_realign.out)
		hsmetrics_run(ABRA2_realign.out)
		mutect2_run(ABRA2_realign.out)
		vardict(ABRA2_realign.out)
		//lofreq(ABRA2_realign.out)
		varscan(ABRA2_realign.out)
		//strelka(ABRA2_realign.out)
		mips_var_caller(ABRA2_realign.out)
		somaticSeq_run(coverage_mosdepth.out.join(mips_var_caller.out.join(mutect2_run.out.join(vardict.out.join(varscan.out.join(ABRA2_realign.out))))))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
