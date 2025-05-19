#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
BED file: ${params.bedfile}.bed
Sequences in:${params.sequences}
"""

include { trimming; coverage_bedtools; coverage_bedtools_uncoll; pair_assembly_pear; pair_assembly_pear as PEAR_UNCOLL; mapping_reads as MAP_UNCOLL; MapBam; GroupReadsByUmi_amplicon } from './mips_mrd.nf'
include { CallMolecularConsensusReadsAmplicon; FilterConsBam; ABRA2_realign as ABRA; hsmetrics_run; mutect2_run; vardict; varscan } from './mips_mrd.nf'
include { hsmetrics_run_uncoll; mutect2_run_uncoll; vardict_uncoll; varscan_uncoll; somaticSeq_run_uncoll } from './mips_mrd.nf'
include { FASTP_UMI } from './scripts/fastp.nf'
include { CUTADAPT_FASTQ; FastqToBam_ASSEMBLED; BAMTOFASTQ } from './scripts/cutadapt.nf'
include { mapping_reads; sam_conversion; GencoreConsensus; ABRA2_realign} from './scripts/mips_gencore.nf'
include { somaticSeq_run; ERROR_CAL } from './scripts/Error_correctn/error.nf'

workflow U2AF1MRD_GENCORE {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.set { samples_ch }

	main:
		CUTADAPT_FASTQ(samples_ch)
		FASTP_UMI(CUTADAPT_FASTQ.out)
		mapping_reads(FASTP_UMI.out)
		sam_conversion(mapping_reads.out)
		GencoreConsensus(sam_conversion.out)
		ABRA2_realign(GencoreConsensus.out)
		coverage_bedtools(ABRA2_realign.out)
		//coverage_bedtools_uncoll(sam_conversion.out)

		// hsmetrics_run(ABRA2_realign.out)
		// hsmetrics_run_uncoll(sam_conversion.out)

		// mutect2_run(ABRA2_realign.out)
		// vardict(ABRA2_realign.out)
		// varscan(ABRA2_realign.out)

		// mutect2_run_uncoll(sam_conversion.out)
		// vardict_uncoll(sam_conversion.out)
		// varscan_uncoll(sam_conversion.out)

		// somaticSeq_run(coverage_bedtools.out.join(mutect2_run.out.join(vardict.out.join(varscan.out.join(ABRA2_realign.out)))))
		// somaticSeq_run_uncoll(coverage_bedtools_uncoll.out.join(mutect2_run_uncoll.out.join(vardict_uncoll.out.join(varscan_uncoll.out.join(sam_conversion.out)))))
}

workflow U2AF1MRD_FGBIO {
	Channel
		.fromPath(params.input)
		.splitCsv(header:false)
		.flatten()
		.set { samples_ch }

	main:
		trimming(samples_ch)
		CUTADAPT_FASTQ(samples_ch)
		FASTP_UMI(CUTADAPT_FASTQ.out)
		pair_assembly_pear(FASTP_UMI.out)
		FastqToBam_ASSEMBLED(pair_assembly_pear.out)
		MapBam(FastqToBam_ASSEMBLED.out)
		GroupReadsByUmi_amplicon(MapBam.out) 
		CallMolecularConsensusReadsAmplicon(GroupReadsByUmi_amplicon.out.filtered_files_ch)
		FilterConsBam(CallMolecularConsensusReadsAmplicon.out) 
		ABRA(FilterConsBam.out)

		PEAR_UNCOLL(trimming.out) | MAP_UNCOLL | sam_conversion

		coverage_bedtools(ABRA.out)
		coverage_bedtools_uncoll(sam_conversion.out)

		hsmetrics_run(ABRA.out)
		hsmetrics_run_uncoll(sam_conversion.out)

		mutect2_run(ABRA.out)
		vardict(ABRA.out)
		varscan(ABRA.out)

		mutect2_run_uncoll(sam_conversion.out)
		vardict_uncoll(sam_conversion.out)
		varscan_uncoll(sam_conversion.out)

		somaticSeq_run(coverage_bedtools.out.join(mutect2_run.out.join(vardict.out.join(varscan.out.join(ABRA.out)))))
		somaticSeq_run_uncoll(coverage_bedtools_uncoll.out.join(mutect2_run_uncoll.out.join(vardict_uncoll.out.join(varscan_uncoll.out.join(sam_conversion.out)))))

		error_ch = somaticSeq_run.out.correctd_excels.collect{ samp_name, file_name -> file_name }
		ERROR_CAL(error_ch)
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the 'Final_Output' directory \n" : "Oops .. something went wrong" )
}