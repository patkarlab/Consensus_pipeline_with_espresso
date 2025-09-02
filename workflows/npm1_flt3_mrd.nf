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
bedfile = file("${params.bedfile}", checkIfExists: true)
illumina_adapters = file("${params.adaptors}", checkIfExists: true)
nextera_adapters = file("${params.nextera_adapters}", checkIfExists: true)
genome_loc = file("${params.genome}", checkIfExists: true)
filt3r_reference = file("${params.filt3r_ref}", checkIfExists: true)
minimap_getitd_reference = file("${params.genome_minimap_getitd}", checkIfExists: true)

// Adapter Trimming and alignment
include { TRIM; MAP } from '../modules/make_bam/bam_npm1_flt3.nf'

// FLT3 ITD detection
include { FILT3R; MINIMAP_GETITD; GETITD; COMBINE_FILT3R_GETITD } from '../modules/flt3_itd/itd_npm1_flt3.nf'

// COVERAGE calculation
include { COVERAGE } from '../modules/bedtools/coverage_npm1_flt3.nf'

// Variant calling
include { VARSCAN } from '../modules/variant_calls/varscan.nf'

// Variant annotation
include { ANNOVAR } from '../modules/annovar/annovar.nf'

// Format output
include { COMBINE_CALLERS } from '../modules/combine_callers/combine_callers.nf'


workflow NPM1_FLT3_MRD {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			def r1 = file("${params.sequences}/${sample}_S*_R1_*.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_S*_R2_*.fastq.gz", checkIfExists: false)

			if (!r1 && !r2) {
				r1 = file("${params.sequences}/${sample}*_R1.fastq.gz", checkIfExists: false)
				r2 = file("${params.sequences}/${sample}*_R2.fastq.gz", checkIfExists: false)
			}
			tuple(sample, r1, r2)
		}
		.set { bam_ch }

	main:
		TRIM(bam_ch, illumina_adapters, nextera_adapters)
		MAP(TRIM.out, genome_loc)
		FILT3R(TRIM.out, filt3r_reference)
		MINIMAP_GETITD(bam_ch, minimap_getitd_reference)
		COVERAGE(MAP.out, bedfile)
		VARSCAN(MAP.out, genome_loc, bedfile)
		ANNOVAR(VARSCAN.out)
		COMBINE_CALLERS(ANNOVAR.out.join(FILT3R.out.join(COVERAGE.out)))
		//COMBINE_FILT3R_GETITD(MINIMAP_GETITD.out.join(COMBINE_CALLERS.out))
}

workflow.onComplete {
	log.info ( workflow.success ? "\n\nDone! Output in the ${params.outdir} directory \n" : "Oops .. something went wrong" )
}
