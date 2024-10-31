#!/usr/bin/env nextflow
nextflow.enable.dsl=2
log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

process SyntheticFastq {
	input:
		val (Sample)
	output:
		tuple val (Sample), file ("*.synreads.bam"), file ("*.synreads.bam.bai")
	script:
	"""
	# Making a synthetic fastq and bam based on the number of reads in the sample bam file
	${params.PosControlScript} ${params.sequences}/${Sample}.bam ${Sample}_inreads
	${params.fastq_bam} ${Sample}_inreads		# This script will output a file named *_inreads.fxd_sorted.bam

	# Merging the Sample bam file with positive control bam file
	cp ${params.sequences}/${Sample}.bam ${Sample}.cons.filtered_merge.bam
	${params.samtools} merge -f ${Sample}.synreads.bam ${Sample}.cons.filtered_merge.bam ${Sample}_inreads.fxd_sorted.bam
	${params.samtools} index ${Sample}.synreads.bam > ${Sample}.synreads.bam.bai
	"""
}

process ABRA2_realign {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.ErrorCorrectd.bam*'
	input:
		tuple val (Sample), file (filt_bam), file (filt_bai)
	output:
		tuple val (Sample), file ("*.ErrorCorrectd.bam"), file ("*.ErrorCorrectd.bam.bai")
	script:
	"""
	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${filt_bam} --out ${Sample}.ErrorCorrectd.bam --ref ${params.genome} --threads 8 --targets ${params.bedfile}.bed --tmpdir ./ > abra.log
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

process mutect2_run {
	maxForks 10
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.mutect2.vcf")
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T MuTect2 -R ${params.genome} -I:tumor ${abra_bam} -o ${Sample}.mutect2.vcf --dbsnp ${params.site2} -L ${params.bedfile}.bed -nct 25 -mbq 20 --allow_potentially_misencoded_quality_scores
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

	sed -i 's/##INFO=<ID=MVD,Number=3,Type=Integer,Description="Calling decision of the 3 algorithms: MuTect, VarScan2, VarDict">/##INFO=<ID=MVD,Number=3,Type=String,Description="Calling decision of the 3 algorithms: MuTect, VarScan2, VarDict">/g' ${Sample}.ErrorCorrectd.vcf

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

process trimming { 
	input:
		val (Sample)
	output:
		//tuple val (Sample), file("${Sample}.R1.trimmed.fastq"), file("${Sample}.R2.trimmed.fastq")
		tuple val (Sample), file("*1P.fq.gz"), file("*2P.fq.gz")
	script:
	"""	
	#/home/programs/ea-utils/clipper/fastq-mcf -o ${Sample}.R1.trimmed.fastq -o ${Sample}.R2.trimmed.fastq -l 53 -k 0 -q 0 /home/diagnostics/pipelines/smMIPS_pipeline/code/functions/preprocess_reads_miseq/smmip_adaptors.fa ${params.sequences}/${Sample}_S*_R1_*.fastq.gz  ${params.sequences}/${Sample}_S*_R2_*.fastq.gz
	trimmomatic PE \
	${params.sequences}/${Sample}_*R1_*.fq.gz ${params.sequences}/${Sample}_*R2_*.fq.gz \
	-baseout ${Sample}.fq.gz \
	ILLUMINACLIP:${params.adaptors}:2:30:10:2:keepBothReads \
	ILLUMINACLIP:${params.nextera_adapters}:2:30:10:2:keepBothReads \
	LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40	
	sleep 5s
	"""
}

process fastqc{
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*'
	input:
		val Sample
	output:
		path "*"
	"""
	${params.fastqc} -o ./ -f fastq ${params.sequences}/${Sample}_*R1_*.fastq.gz ${params.sequences}/${Sample}_*R2_*.fastq.gz
	"""
}

process pair_assembly_pear {
	input:
		tuple val (Sample), file(paired_forward), file(paired_reverse)
	output:
		tuple val (Sample), file("*.assembled.fastq")
	script:
	"""
	${params.pear_path} -f ${paired_forward} -r ${paired_reverse} -o ${Sample} -n 53 -j 25
	"""
}

process mapping_reads {
	input:
		tuple val (Sample), file (pairAssembled)
	output:
		tuple val (Sample), file ("*.sam")
	script:
	"""
	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200" -M -t 20 ${params.genome} ${pairAssembled} > ${Sample}.sam
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

process coverage_mosdepth_uncoll {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.regions.bed'
	input:
		tuple val (Sample), file (abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}_uncollaps.regions.bed"), file ("${Sample}_exon_uncollaps.regions.bed")
	script:
	"""
	${params.mosdepth_script} ${abra_bam} ${Sample}_uncollaps ${params.bedfile}.bed
	${params.mosdepth_script} ${abra_bam} ${Sample}_exon_uncollaps ${params.bedfile2}.bed
	"""
}

process hsmetrics_run_uncoll {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_uncoll_hsmetrics.txt'
	input:
		tuple val(Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*_uncoll_hsmetrics.txt")
	script:
	"""
	${params.java_path}/java -jar ${params.picard_path} CollectHsMetrics I= ${abra_bam} O= ${Sample}_uncoll_hsmetrics.txt BAIT_INTERVALS= ${params.bedfile}.interval_list TARGET_INTERVALS= ${params.bedfile}.interval_list R= ${params.genome} VALIDATION_STRINGENCY=LENIENT
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

process mutect2_run_uncoll {
	maxForks 10
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("${Sample}.mutect2.vcf")
	script:
	"""
	${params.java_path}/java -Xmx10G -jar ${params.GATK42_path} Mutect2 -R ${params.genome} -I:tumor ${abra_bam} -O ${Sample}.mutect2.vcf -L ${params.bedfile}.bed -mbq 20
	${params.java_path}/java -Xmx10G -jar ${params.GATK42_path} FilterMutectCalls -V ${Sample}.mutect2.vcf --stats ${Sample}.mutect2.vcf.stats -O ${Sample}_filtered.vcf -R ${params.genome} --unique-alt-read-count 3 --min-median-base-quality 20 --min-median-mapping-quality 30
	#perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.mutect2.vcf --outfile ${Sample}.mutect2.avinput -allsample -withfreq --includeinfo
	#perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.mutect2.avinput --out ${Sample}.mutect2.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	#if [ -s ${Sample}.mutect2.somaticseq.hg19_multianno.csv ]; then
	#	python3 ${PWD}/scripts/somaticseqoutput-format_v2_mutect2.py ${Sample}.mutect2.somaticseq.hg19_multianno.csv mutect2_.csv
	#else
	#	touch mutect2_.csv		
	#fi	
	"""
}

process vardict_uncoll {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.0001 -r 8 -N ${Sample} -b ${abra_bam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf
	"""
}

process varscan_uncoll {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -d 1000000 -A -a -l ${params.bedfile}.bed -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup -min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup -min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
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
		SyntheticFastq(samples_ch)
		ABRA2_realign(SyntheticFastq.out)
		coverage_mosdepth(ABRA2_realign.out)
		hsmetrics_run(ABRA2_realign.out)
		CNS_filegen(ABRA2_realign.out)
		espresso(CNS_filegen.out.collect())
		mutect2_run(ABRA2_realign.out)
		vardict(ABRA2_realign.out)
		//lofreq(ABRA2_realign.out)
		varscan(ABRA2_realign.out)
		//strelka(ABRA2_realign.out)
		somaticSeq_run(coverage_mosdepth.out.join(mutect2_run.out.join(vardict.out.join(varscan.out.join(ABRA2_realign.out)))))

		//Uncollapsed data analysis
		trimming(samples_ch)
		fastqc(samples_ch)
		pair_assembly_pear(trimming.out) | mapping_reads | sam_conversion
		coverage_mosdepth_uncoll(sam_conversion.out)
		hsmetrics_run_uncoll(sam_conversion.out)
		minimap_getitd(sam_conversion.out)
		mutect2_run_uncoll(sam_conversion.out)
		vardict_uncoll(sam_conversion.out)
		//lofreq_uncoll(MapBam.out)
		varscan_uncoll(sam_conversion.out)
		somaticSeq_run_uncoll(coverage_mosdepth_uncoll.out.join(mutect2_run_uncoll.out.join(vardict_uncoll.out.join(varscan_uncoll.out.join(sam_conversion.out)))))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
