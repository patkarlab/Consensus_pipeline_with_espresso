#!/usr/bin/env nextflow
nextflow.enable.dsl=2

log.info """
STARTING PIPELINE
=*=*=*=*=*=*=*=*=
Sample list: ${params.input}
Sequences in:${params.sequences}
"""

// include { ADD_UMI } from './modules/fastq_mcf/preprocess/main.nf'



workflow MRD_PROBE {
	Channel.fromPath(params.input)
		.splitCsv(header:false)
		.map { row ->
			def sample = row[0].trim()
			def r1 = file("${params.sequences}/${sample}_R1.fastq.gz", checkIfExists: false)
			def r2 = file("${params.sequences}/${sample}_R2.fastq.gz", checkIfExists: false)
			def umi = file("${params.sequences}/${sample}_UMI_R1.fastq.gz", checkIfExists: false)
			tuple(sample_short, r1, r2, umi)
		}
		.view()
		// .set { bam_ch }

	// main:
	// ADD_UMI(bam_ch)
}