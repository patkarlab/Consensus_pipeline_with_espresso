#!/usr/bin/env nextflow
nextflow.enable.dsl=2


genome_loc = file("${params.genome}", checkIfExists: true )
genome_dir = file("${genome_loc.parent}", checkIfExists: true)
genome_fasta = file("${genome_loc.name}")
bedfile = file("${params.bedfile}.bed", checkIfExists: true)
interval_list = file("${params.bedfile}.interval_list", checkIfExists: true)
truseq_adapters = file("${params.adaptors}", checkIfExists: true)
nextera_adapters = file("${params.nextera_adapters}", checkIfExists: true)


include { COVERAGE as COVERAGE_COLL; COVERAGE as  COVERAGE_UNCOLL } from './modules/cebpa_amplicon/main.nf'
include { HSMETRICS as HSMETRICS_COLL; HSMETRICS as HSMETRICS_UNCOLL } from './modules/cebpa_amplicon/main.nf'
include { MUTECT2 as MUTECT2_COLL; MUTECT2 as MUTECT2_UNCOLL} from './modules/cebpa_amplicon/main.nf'
include { VARDICT as VARDICT_COLL; VARDICT as VARDICT_UNCOLL} from './modules/cebpa_amplicon/main.nf'
include { VARSCAN as VARSCAN_COLL; VARSCAN as VARSCAN_UNCOLL} from './modules/cebpa_amplicon/main.nf'
include { ANNOVAR as ANNOVAR_COLL_MUTECT2; ANNOVAR as ANNOVAR_COLL_VARDICT; ANNOVAR as ANNOVAR_COLL_VARSCAN } from './modules/cebpa_amplicon/main.nf'
include { ANNOVAR as ANNOVAR_UNCOLL_MUTECT2; ANNOVAR as ANNOVAR_UNCOLL_VARDICT; ANNOVAR as ANNOVAR_UNCOLL_VARSCAN } from './modules/cebpa_amplicon/main.nf'
include { COMBINE_CALLERS as COMBINE_CALLERS_COLL; COMBINE_CALLERS as COMBINE_CALLERS_UNCOLL} from './modules/cebpa_amplicon/main.nf'


workflow CEBPA_MRD {
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
		.set { samples_ch }

	main:
		TRIM(samples_ch, truseq_adapters, nextera_adapters)

		FASTQTOBAM(trimming.out) 
		MAPBAM(FASTQTOBAM.out, genome_dir, genome_fasta) 
		GROUPREADSBYUMI(MAPBAM.out) 
		CALLMOLCONSREADS(GROUPREADSBYUMI.out.grouped_bam_ch, genome_dir, genome_fasta) 
		FILTERCONSBAM(CALLMOLCONSREADS.out, genome_dir, genome_fasta)
		ABRA2_REALIGN(FILTERCONSBAM.out, bedfile, genome_dir, genome_fasta)

		MAPBAM_UNCOLL(TRIM.out)

		COVERAGE_COLL(ABRA2_REALIGN.out, genome_dir, genome_fasta, bedfile, collapsed)
		COVERAGE_UNCOLL(MAPBAM_UNCOLL.out, genome_dir, genome_fasta, bedfile, uncollapsed)

		HSMETRICS_COLL(ABRA2_REALIGN.out, genome_dir, genome_fasta, interval_list, collapsed)
		HSMETRICS_UNCOLL(MAPBAM_UNCOLL.out, genome_dir, genome_fasta, interval_list, uncollapsed)

		MUTECT2_COLL(ABRA2_REALIGN.out, genome_dir, genome_fasta, bedfile)
		VARDICT_COLL(ABRA2_REALIGN.out, genome_dir, genome_fasta, bedfile)
		VARSCAN_COLL(ABRA2_REALIGN.out, genome_dir, genome_fasta, bedfile)

		//MUTECT2_UNCOLL(MAPBAM_UNCOLL.out)
		//VARDICT_UNCOLL(MAPBAM_UNCOLL.out)
		//VARSCAN_UNCOLL(MAPBAM_UNCOLL.out)    

		ANNOVAR_COLL_MUTECT2(MUTECT2_COLL.out, mutect2)
		ANNOVAR_COLL_VARDICT(VARDICT_COLL.out, vardict)
		ANNOVAR_COLL_VARSCAN(VARSCAN_COLL.out, varscan)

		//ANNOVAR_UNCOLL_MUTECT2(MUTECT2_UNCOLL.out, mutect2)
		//ANNOVAR_UNCOLL_VARDICT(VARDICT_UNCOLL.out, vardict)
		//ANNOVAR_UNCOLL_VARSCAN(VARSCAN_UNCOLL.out, varscan)

		COMBINE_CALLERS_COLL(ANNOVAR_COLL_MUTECT2.out.join(ANNOVAR_COLL_VARDICT.out.join(ANNOVAR_COLL_VARSCAN.out.join(COVERAGE_COLL.out))), collapsed)
		//COMBINE_CALLERS_UNCOLL(ANNOVAR_UNCOLL_MUTECT2.out.join(ANNOVAR_UNCOLL_VARDICT.out.join(ANNOVAR_UNCOLL_VARSCAN.out.join(COVERAGE_UNCOLL.out))), uncollapsed)
		
}