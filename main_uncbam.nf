#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""
process FastqToBam { 
	input:
		val Sample
	output:
		tuple val (Sample), file("*.unmapped.bam")
	script:
	"""	
	java -Xmx12g -jar ${params.fgbio_path} --compression 1 --async-io FastqToBam --input ${params.sequences}/${Sample}_*_R1_*.fastq.gz ${params.sequences}/${Sample}_*_R2_*.fastq.gz ${params.sequences}/${Sample}_*_R3_*.fastq.gz --read-structures +T +M +T --sample ${Sample} --library ${Sample} --output ${Sample}.unmapped.bam
	sleep 5s
	"""
}

process MapBam {	
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.sorted.bam*'
	input:
		tuple val (Sample), file(unmapped_bam)
	output:
		tuple val (Sample), file("*.sorted.bam"), file ("*.sorted.bam.bai")
	script:
	"""
	${params.samtools} fastq ${unmapped_bam} | bwa mem -t 16 -p -K 150000000 -Y ${params.genome} - | java -Xmx4g -jar ${params.fgbio_path} --compression 1 --async-io ZipperBams --unmapped ${unmapped_bam} --ref ${params.genome} --output ${Sample}.mapped.bam

	${params.samtools} sort ${Sample}.mapped.bam -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam > ${Sample}.sorted.bam.bai
	"""
}

process GroupReadsByUmi {
	input:
		tuple val (Sample), file(mapped_bam), file(mapped_bai)
	output: 
		tuple val (Sample), file ("*.grouped.bam")	
	script:
	"""
	java -Xmx8g -jar ${params.fgbio_path} --compression 1 --async-io GroupReadsByUmi --input ${mapped_bam} --strategy Adjacency --edits 1 --output ${Sample}.grouped.bam --family-size-histogram ${Sample}.tag-family-sizes.txt
	"""
}

process CallMolecularConsensusReads {
	input:
		tuple val (Sample), file (grouped_bam)
	output:
		tuple val (Sample), file ("*.cons.unmapped.bam")
	script:
	"""
	java -Xmx4g -jar ${params.fgbio_path} --compression 0 CallMolecularConsensusReads --input ${grouped_bam} --output /dev/stdout --min-reads 3 --min-input-base-quality 20 --threads 4 | java -Xmx4g -jar ${params.fgbio_path} --compression 1 FilterConsensusReads --input /dev/stdin --output ${Sample}.cons.unmapped.bam --ref ${params.genome} --min-reads 3 --min-base-quality 45 --max-base-error-rate 0.2
	"""
}

process FilterConsBam {
	input:
		tuple val (Sample), file (cons_unmapped_bam)
	output:
		tuple val (Sample), file ("*.cons.filtered.bam")
	script:
	"""
	${params.samtools} fastq ${cons_unmapped_bam} | bwa mem -t 16 -p -K 150000000 -Y ${params.genome} - | java -Xmx12g -jar ${params.fgbio_path} --compression 0 --async-io ZipperBams --unmapped ${cons_unmapped_bam} --ref ${params.genome} --tags-to-reverse Consensus --tags-to-revcomp Consensus | ${params.samtools} sort --threads 8 -o ${Sample}.cons.filtered.bam
	"""
}

process SyntheticFastq {
	input:
		tuple val (Sample), file (cons_filt_bam)
	output:
		tuple val (Sample), file ("*.synreads.bam"), file ("*.synreads.bam.bai")
	script:
	"""
	# Making a synthetic fastq and bam based on the number of reads in the sample bam file
	${params.PosControlScript} ${cons_filt_bam} ${Sample}_inreads
	${params.fastq_bam} ${Sample}_inreads		# This script will output a file named *_inreads.fxd_sorted.bam

	# Merging the Sample bam file with positive control bam file
	mv ${cons_filt_bam} ${Sample}.cons.filtered_merge.bam
	${params.samtools} merge -f ${Sample}.synreads.bam ${Sample}.cons.filtered_merge.bam ${Sample}_inreads.fxd_sorted.bam
	${params.samtools} index ${Sample}.synreads.bam > ${Sample}.synreads.bam.bai
	"""
}

process ABRA2_realign {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.abra.bam*'
	input:
		tuple val (Sample), file (filt_bam), file (filt_bai)
	output:
		tuple val (Sample), file ("*.abra.bam"), file ("*.abra.bam.bai")
	script:
	"""
	${params.java_path}/java -Xmx16G -jar ${params.abra2_path}/abra2-2.23.jar --in ${filt_bam} --out ${Sample}.abra.bam --ref ${params.genome} --threads 8 --targets ${params.bedfile}.bed --tmpdir ./ > abra.log
	${params.samtools} index ${Sample}.abra.bam > ${Sample}.abra.bam.bai
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
	${params.samtools} mpileup -s -x -BQ 0 -q 1 -d 1000000 --skip-indels -f ${params.genome} ${abra_bam} > ${Sample}.mpileup

	# These parmeters were taken from the Espresso paper(Abelson et al., Sci. Adv. 2020; 6 : eabe3722)
	# A stringent (lower) p value removes all the snps hence, p-value is 1
	${params.java_path}/java -jar ${params.varscan_path} pileup2cns ${Sample}".mpileup" --variants SNP --min-coverage 10 --min-reads2 1 --min-avg-qual 30 --min-var-freq 0.0001 --p-value 1 --strand-filter 0 > ${Sample}".cns"

	grep -v '^chrM' ${Sample}".cns" > ${Sample}".nochrM.cns"
	${params.filter_cns} ${Sample}".nochrM.cns" ${Sample}".final.cns"
	"""
}

process coverage_mosdepth {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*cov_uncbam.regions.bed'
	input:
		tuple val (Sample), file (abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	${params.mosdepth_script} ${abra_bam} ${Sample}_cov ${params.bedfile}.bed
	mv ${Sample}_cov.regions.bed ${Sample}_cov_uncbam.regions.bed
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
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*_getitd_uncbam'
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		path "*_getitd_uncbam"
	script:
	"""
	${params.samtools} sort ${abra_bam} -o ${Sample}.sorted.bam
	${params.samtools} index ${Sample}.sorted.bam
	${params.samtools} view ${Sample}.sorted.bam -b -h chr13 > ${Sample}.chr13.bam

	samtools sort ${Sample}.chr13.bam -o temp.bam 
	samtools index temp.bam > temp.bam.bai

	${params.bedtools} bamtofastq -i temp.bam -fq ${Sample}_chr13.fastq
	python ${params.get_itd_path}/getitd.py -reference ${params.get_itd_path}/anno/amplicon.txt -anno ${params.get_itd_path}/anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 ${Sample} ${Sample}_chr13.fastq
	mv ${Sample}_getitd ${Sample}_getitd_uncbam
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
	#${params.java_path}/java -Xmx10G -jar ${params.GATK38_path} -T MuTect2 -R ${params.genome} -I:tumor ${abra_bam} -o ${Sample}.mutect2.vcf --dbsnp ${params.site2} -L ${params.bedfile}.bed -nct 25 -contamination 0.02 -mbq 30
	${params.samtools} view -bs 40.1 ${abra_bam} > subsampled_01.bam
	${params.samtools} index subsampled_01.bam
	${params.mutect2} ${params.java_path} ${params.GATK38_path} ${params.genome} subsampled_01.bam ${Sample}.mutect2.vcf ${params.site2} ${params.bedfile}.bed
	"""
}

process vardict {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.vardict.vcf")
	script:
	"""
	VarDict -G ${params.genome} -f 0.03 -N ${Sample} -b ${abra_bam} -c 1 -S 2 -E 3 -g 4 ${params.bedfile}.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.03 > ${Sample}.vardict.vcf
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
	${params.lofreq_path} call -b dynamic -C 50 -a 0.00005 -q 30 -Q 30 -m 50 -f ${params.genome} -l ${params.bedfile}.bed -o ${Sample}.lofreq.vcf ${Sample}.lofreq.bam
	${params.lofreq_path} filter -a 0.005 -i ${Sample}.lofreq.vcf -o ${Sample}.lofreq.filtered.vcf
	"""
}

process varscan {
	input:
		tuple val (Sample), file(abra_bam), file (abra_bai)
	output:
		tuple val (Sample), file ("*.varscan.vcf")
	script:
	"""
	${params.samtools} mpileup -f ${params.genome} ${abra_bam} > ${Sample}.mpileup
	${params.java_path}/java -jar ${params.varscan_path} mpileup2snp ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_snp.vcf
	${params.java_path}/java -jar ${params.varscan_path} mpileup2indel ${Sample}.mpileup --min-coverage 10 --min-reads2 5 --min-avg-qual 15 --min-var-freq 0.003 --p-value 1e-4 --output-vcf 1 > ${Sample}.varscan_indel.vcf
	bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
	bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_snp.vcf.gz
	${params.bcftools_path} index -t ${Sample}.varscan_indel.vcf.gz
	${params.bcftools_path} concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
	"""
}

process somaticSeq_run {
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.somaticseq_uncbam.vcf'
	publishDir "$PWD/Final_Output/${Sample}/", mode: 'copy', pattern: '*.somaticseq_uncbam.csv'
	input:
		tuple val (Sample), file (lofreqVcf), file (vardictVcf), file (varscanVcf), file(finalBam), file (finalBamBai)
	output:
		tuple val (Sample), file ("*")
	script:
	"""
	samtools sort ${finalBam} -o temp.bam
	samtools index temp.bam > temp.bam.bai

	somaticseq_parallel.py --output-directory ${Sample}.somaticseq --genome-reference ${params.genome} --inclusion-region ${params.bedfile}.bed --threads 25 --algorithm xgboost  --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file temp.bam --vardict-vcf ${vardictVcf} --varscan-vcf ${varscanVcf} --lofreq-vcf ${lofreqVcf} --sample-name ${Sample}

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sSNV.vcf ${Sample}.somaticseq/somaticseq_snv.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_snv.vcf > ${Sample}.somaticseq/somaticseq_snv.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_snv.vcf.gz

	${params.vcf_sorter_path} ${Sample}.somaticseq/Consensus.sINDEL.vcf ${Sample}.somaticseq/somaticseq_indel.vcf
	bgzip -c ${Sample}.somaticseq/somaticseq_indel.vcf > ${Sample}.somaticseq/somaticseq_indel.vcf.gz
	${params.bcftools_path} index -t ${Sample}.somaticseq/somaticseq_indel.vcf.gz

	${params.bcftools_path} concat -a ${Sample}.somaticseq/somaticseq_snv.vcf.gz ${Sample}.somaticseq/somaticseq_indel.vcf.gz -o ${Sample}.somaticseq.vcf

	perl ${params.annovarLatest_path}/convert2annovar.pl -format vcf4 ${Sample}.somaticseq.vcf  --outfile ${Sample}.somaticseq.avinput --withzyg --includeinfo
	perl ${params.annovarLatest_path}/table_annovar.pl ${Sample}.somaticseq.avinput --out ${Sample}.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 ${params.annovarLatest_path}/humandb/ --xreffile ${params.annovarLatest_path}/example/gene_fullxref.txt

	if [ -s ${Sample}.somaticseq.hg19_multianno.csv ]; then
		python3 ${params.format_somaticseq_script} ${Sample}.somaticseq.hg19_multianno.csv ${Sample}.somaticseq.csv
	else
		touch ${Sample}.somaticseq.csv
	fi

	mv ${Sample}.somaticseq.vcf ${Sample}.somaticseq_uncbam.vcf
	mv ${Sample}.somaticseq.csv ${Sample}.somaticseq_uncbam.csv
	"""
}

workflow MRD {

	Channel
		.fromPath(params.input)

		.splitCsv(header:false)
		.flatten()
		.map{ it }
		.set { samples_ch }

	main:
		FastqToBam(samples_ch)
		MapBam(FastqToBam.out) 		
		//GroupReadsByUmi(MapBam.out) 
		//CallMolecularConsensusReads(GroupReadsByUmi.out) 
		//FilterConsBam(CallMolecularConsensusReads.out) 
		//SyntheticFastq(FilterConsBam.out)
		//ABRA2_realign(SyntheticFastq.out)
		//CNS_filegen(ABRA2_realign.out)
		//coverage_mosdepth(ABRA2_realign.out)
		//coverage_mosdepth(MapBam.out)
		//hsmetrics_run(ABRA2_realign.out)
		//minimap_getitd(MapBam.out)
		//vardict(ABRA2_realign.out)
		//lofreq(ABRA2_realign.out)
		//varscan(ABRA2_realign.out)
		vardict(MapBam.out)
		lofreq(MapBam.out)
		varscan(MapBam.out)
		//somaticSeq_run(mutect2_run.out.join(varscan.out.join(ABRA2_realign.out)))
		//somaticSeq_run(mutect2_run.out.join(varscan.out.join(MapBam.out)))
		somaticSeq_run(lofreq.out.join(vardict.out.join(varscan.out.join(MapBam.out))))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}
