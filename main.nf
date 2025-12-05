#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""
// file paths
bed_file = file("${params.bedfile}.bed", checkIfExists: true )
genome_file = file("${params.genome}", checkIfExists: true)
matrix_file = file("${params.matrix}", checkIfExists: true)

include { FASTQTOBAM; TRIM } from './modules/make_bam/bam'
include { FastqToBam; MapBam; GroupReadsByUmi; CallMolecularConsensusReads } from './modules/fgbio/collapsed_bam.nf'
include { COMBINEBAMS; FilterConsBam; SyntheticFastq } from './modules/fgbio/collapsed_bam.nf'
include { ABRA2_realign } from './modules/abra/abra.nf'
include { CNS_FILEGEN; ESPRESSO } from './modules/espresso/espresso_error.nf'
include { COVERAGE_BEDTOOLS; COVERAGE_BEDTOOLS_UNCOLL } from './modules/bedtools/coverage.nf'
include { HSMETRICS; HSMETRICS_UNCOLL } from './modules/picard/hsmetrics.nf'
include { GETITD } from './modules/flt3_itd/get_itd.nf'
include { VARDICT as VARDICT_COLLAPSE; VARDICT as VARDICT_UNCOLLAPSE } from './modules/variant_calls/variant_call.nf'
include { VARSCAN as VARSCAN_COLLAPSE; VARSCAN as VARSCAN_UNCOLLAPSE } from './modules/variant_calls/variant_call.nf'
include { MUTECT2 as MUTECT2_COLLAPSE; MUTECT2 as MUTECT2_UNCOLLAPSE } from './modules/variant_calls/variant_call.nf'
include { DEEPSOMATIC as DEEPSOMATIC_COLLAPSE; DEEPSOMATIC as DEEPSOMATIC_UNCOLLAPSE } from './modules/variant_calls/variant_call.nf'
include { somaticSeq_run; somaticSeq_run_uncoll } from './modules/somaticseq/somaticseq.nf'
include { ERROR_CAL; ERROR_MODEL } from './modules/Error_correctn/error.nf'
include { SPLITBAM } from './modules/fgbio/collapsed_bam.nf'
include { ANNOVAR_UNCOLLAPSE; ANNOVAR } from './modules/variant_calls/variant_call.nf'

workflow MIPS_MRD {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			def sample_short = (sample =~ /_(.*?)-Pool/)[0][1]
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)
			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample_short, r1, r2)
		}
		.set { bam_ch }
		// .view()
	main:
	 	// Uncollapsed data analysis
		FASTQTOBAM(bam_ch)
		COVERAGE_BEDTOOLS_UNCOLL(FASTQTOBAM.out.final_bams_ch)
		GETITD (FASTQTOBAM.out.final_bams_ch)
		VARDICT_UNCOLLAPSE(FASTQTOBAM.out.final_bams_ch)
		VARSCAN_UNCOLLAPSE(FASTQTOBAM.out.final_bams_ch)
		MUTECT2_UNCOLLAPSE(FASTQTOBAM.out.final_bams_ch)
		ANNOVAR_UNCOLLAPSE(COVERAGE_BEDTOOLS_UNCOLL.out.join(MUTECT2_UNCOLLAPSE.out.join(VARDICT_UNCOLLAPSE.out.join(VARSCAN_UNCOLLAPSE.out))))

		// Collapsed data analysis
		// FastqToBam(FASTQTOBAM.out.trim_ch)
		// MapBam(FastqToBam.out)
		// SPLITBAM(MapBam.out, bed_file)
		// GroupReadsByUmi(SPLITBAM.out.chrwise_bamlist) 
		// CallMolecularConsensusReads(GroupReadsByUmi.out.grouped_bam_ch)
		// COMBINEBAMS (CallMolecularConsensusReads.out.consensus_bam_ch) 
		// FilterConsBam(COMBINEBAMS.out.combined_bam_ch) 
		// SyntheticFastq(FilterConsBam.out)
		// ABRA2_realign(SyntheticFastq.out)
		// COVERAGE_BEDTOOLS(ABRA2_realign.out)
		// VARDICT_COLLAPSE(ABRA2_realign.out)
		// VARSCAN_COLLAPSE(ABRA2_realign.out)
		// MUTECT2_COLLAPSE(ABRA2_realign.out)
		// ANNOVAR(COVERAGE_BEDTOOLS.out.join(MUTECT2_COLLAPSE.out.join(VARDICT_COLLAPSE.out.join(VARSCAN_COLLAPSE.out))))

		// error_ch = somaticSeq_run.out.correctd_excels.collect{ samp_name, file_name -> file_name }
		// ERROR_CAL(error_ch)
		// ERROR_MODEL(ANNOVAR.out, matrix_file)
}

workflow BNC_NORM {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			def sample_short = (sample =~ /_(.*?)-Pool/)[0][1]
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)
			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample_short, r1, r2)
		}
		// .set { bam_ch }
		.view()
	// main:
		// TRIM(bam_ch)
		// FastqToBam(TRIM.out)
		// MapBam(FastqToBam.out)
		// SPLITBAM(MapBam.out, bed_file)
		// GroupReadsByUmi(SPLITBAM.out.chrwise_bamlist) 
		// CallMolecularConsensusReads(GroupReadsByUmi.out.grouped_bam_ch)
		// COMBINEBAMS (CallMolecularConsensusReads.out.consensus_bam_ch) 
		// FilterConsBam(COMBINEBAMS.out.combined_bam_ch)
		// SyntheticFastq(FilterConsBam.out)
		// ABRA2_realign(SyntheticFastq.out)
		// MUTECT2_COLLAPSE(ABRA2_realign.out)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
	log.info ("Completed at: ${workflow.complete}")
	log.info ("Total time taken: ${workflow.duration}")
}