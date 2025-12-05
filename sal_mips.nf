#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
Sequences in:${params.sequences}
"""

// file paths
genome_file = file("${params.genome}", checkIfExists: true)
matrix_file = file("${params.matrix}", checkIfExists: true)
adaptors_file = file("${params.smmips_adaptors}", checkIfExists:true)
index_files = file("${params.genome_dir}/${params.ind_files}.*")
txt_file = file("${params.smMIPS_txt_file}", checkIfExists:true)
txtR_file = file("${params.smMIPS_txtR_file}", checkIfExists:true)
bin_script_files = file("${params.bin_scripts}/*", checkIfExists:true)
// getitd_amplicon = file("${params.get_itd_path}/anno/amplicon.txt", checkIfExists:true)
// getitd_amplicon_kayser = file("${params.get_itd_path}/anno/amplicon_kayser.tsv", checkIfExists:true)
sam_forward = params.forward_sam
sam_reverse = params.reverse_sam
annotate_snp = params.snp
annotate_indel = params.indel
beta_matrix = file("${params.matrix}", checkIfExists:true)
snp_list = file("${params.snp_filter_list}", checkIfExists:true)
indel_list = file("${params.indel_filter_list}", checkIfExists:true)

include { PREPROCESS } from './modules/fastq_mcf/preprocess/main.nf'
include { PEAR_SELF_ASSEMBLE } from './modules/pear/assemble/main.nf'
include { RENAME_SMMIPS } from './modules/rename/main.nf'
include { ALIGNMENT } from './modules/bwa/align/main.nf'
include { SPLITSAM } from './modules/samtools/view/main.nf'
include { BAMTOFASTQ } from './modules/bedtools/bamtofastq/main.nf'
// include { GETITD } from './modules/flt3_itd/get_itd/main.nf'
include { SPLIT_MIPS as SPLIT_MIPS_F ; SPLIT_MIPS as SPLIT_MIPS_R } from './modules/split_smmips/main.nf'
include { CALL_CONSENSUS as CALL_CONSENSUS_F; CALL_CONSENSUS as CALL_CONSENSUS_R} from './modules/fgbio/CallMolecularConsensus/main.nf'
include { SAMTOOLS as SAMTOOLS_F; SAMTOOLS as SAMTOOLS_R } from './modules/samtools/sort/main.nf'
include { SAMTOFASTQ as SAMTOFASTQ_F; SAMTOFASTQ as SAMTOFASTQ_R } from './modules/gatk/SamToFastq/main.nf'
include { MPILEUP as MPILEUP_F; MPILEUP as MPILEUP_R} from './modules/bwa/mpileup/main.nf'
include { VARCALL as VARCALL_F; VARCALL as VARCALL_R} from './modules/samtools/variant_calls/main.nf'
include { COMBINE_VCF } from './modules/grep/combine_vcf/main.nf'
include { COMBINE_MIPS_COUNTS } from './modules/cat/combine_mipcount/main.nf'
include { ANNOVAR as ANNOVAR_SNP; ANNOVAR as ANNOVAR_INDEL } from './modules/annovar/annotate/main.nf'
include { FORMAT_ANNOVAR_SNP; FORMAT_ANNOVAR as FORMAT_ANNOVAR_INDEL } from './modules/python/format/main.nf'
include { ERROR_CORRECTN } from './modules/beta_dist_model/error_correctn/main.nf'

workflow MIPS {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			// def sample_short = (sample =~ /_(.*?)-Pool/)[0][1]
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)
			if (!r1 || !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample, r1, r2)
		}
		.set { bam_ch }
	main:
	PREPROCESS(bam_ch, adaptors_file)
	PEAR_SELF_ASSEMBLE(PREPROCESS.out)
	RENAME_SMMIPS(PEAR_SELF_ASSEMBLE.out)
	ALIGNMENT(RENAME_SMMIPS.out, genome_file, index_files)
	SPLITSAM(ALIGNMENT.out.bam_out_ch)
	// BAMTOFASTQ(ALIGNMENT.out.chr13_out_ch)
	// GETITD(BAMTOFASTQ.out, getitd_amplicon, getitd_amplicon_kayser)
	SPLIT_MIPS_F(SPLITSAM.out, txt_file, sam_forward, bin_script_files)
	SPLIT_MIPS_R(SPLITSAM.out, txtR_file, sam_reverse, bin_script_files)
	CALL_CONSENSUS_F(SPLIT_MIPS_F.out.nosingles_sam_ch)
	CALL_CONSENSUS_R(SPLIT_MIPS_R.out.nosingles_sam_ch)
	SAMTOOLS_F(CALL_CONSENSUS_F.out)
	SAMTOOLS_R(CALL_CONSENSUS_R.out)
	SAMTOFASTQ_F(SAMTOOLS_F.out)
	SAMTOFASTQ_R(SAMTOOLS_R.out)
	MPILEUP_F(SAMTOFASTQ_F.out, genome_file, index_files)
	MPILEUP_R(SAMTOFASTQ_R.out, genome_file, index_files)
	VARCALL_F(MPILEUP_F.out.mpileup_ch.join(CALL_CONSENSUS_F.out), txt_file, sam_forward)
	VARCALL_R(MPILEUP_R.out.mpileup_ch.join(CALL_CONSENSUS_R.out), txtR_file, sam_reverse)
	COMBINE_VCF(VARCALL_F.out.join(VARCALL_R.out))
	COMBINE_MIPS_COUNTS(SPLIT_MIPS_F.out.mip_count_ch.join(SPLIT_MIPS_R.out.mip_count_ch))
	ANNOVAR_SNP(COMBINE_VCF.out.vcf_point_ch, annotate_snp)
	ANNOVAR_INDEL(COMBINE_VCF.out.vcf_indel_ch, annotate_indel)
	FORMAT_ANNOVAR_SNP(ANNOVAR_SNP.out)
	FORMAT_ANNOVAR_INDEL(ANNOVAR_INDEL.out, indel_list)
	ERROR_CORRECTN(FORMAT_ANNOVAR_SNP.out, annotate_snp, beta_matrix)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
	log.info ("Completed at: ${workflow.complete}")
	log.info ("Total time taken: ${workflow.duration}")
}