#!/usr/bin/nextflow

process FastqToBam {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file(R1_trim_fastq), file(R2_trim_fastq)
	output:
		tuple val (Sample), file ("*.unmapped.bam")
	script:	
	"""
	java -Xmx${task.memory.toGiga()}g -jar "/home/programs/fgbio/fgbio-2.0.1.jar" --compression 1 --async-io FastqToBam \
	--input ${R1_trim_fastq} ${R2_trim_fastq} --read-structures 4M+T 4M+T \
	--sample ${Sample} --library ${Sample} --output ${Sample}.unmapped.bam
	"""
}

process MapBam {
	tag "${Sample}"
	label 'process_low'	
	input:
		tuple val (Sample), file(unmapped_bam)
	output:
		tuple val (Sample), file("*.mapped.bam")
	script:
	"""
	# Align the data with bwa and recover headers and tags from the unmapped BAM
	${params.samtools} fastq ${unmapped_bam} | bwa mem -t ${task.cpus} -p -K 150000000 -Y ${params.genome} - | \
	java -Xmx${task.memory.toGiga()}g -jar ${params.fgbio_path} --compression 1 --async-io ZipperBams \
	--unmapped ${unmapped_bam} --ref ${params.genome} --output ${Sample}.mapped.bam
	"""
}

process SPLITBAM {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file (mapbam)
		path (bedfile)
	output:
		tuple val (Sample), path ("${Sample}_*_primary.bam"), emit: chrwise_bamlist
	script:
	"""
	samtools sort -@ ${task.cpus} ${mapbam} -o ${Sample}_sortd.bam 
	samtools index ${Sample}_sortd.bam > ${Sample}_sortd.bam.bai
	samtools idxstats ${Sample}_sortd.bam | cut -f 1 | grep -v '*' > chromosomes.txt
	for i in `cat chromosomes.txt`
	do 
		#samtools view -@ ${task.cpus} -b ${Sample}_sortd.bam "\${i}" -o ${Sample}_\${i}_split.bam
		# Obtain the read names (headers) aligned to each chromosome
		samtools view -@ ${task.cpus} ${Sample}_sortd.bam "\${i}" | cut -f1 | sort -u > \${i}_names.txt
		# Mate rescue: extract *all* alignments for those names
		samtools view -@ ${task.cpus} -N \${i}_names.txt -b ${Sample}_sortd.bam > ${Sample}_\${i}_rescued.bam
		# 3) Keep only primaries (drop secondary=256 and supplementary=2048)
		samtools view -@ ${task.cpus} -F 2304 -b ${Sample}_\${i}_rescued.bam > ${Sample}_\${i}_primary.bam
	done

	#counter=0
	#while read chr start stop region; do
	#	samtools view -@ ${task.cpus} -P -b ${Sample}_sortd.bam \${chr}:\${start}-\${stop} > ${Sample}_split_\${counter}.bam
	#	#samtools sort -@ ${task.cpus} ${Sample}_split_\${counter}.bam -o ${Sample}_split_\${counter}_sort.bam
	#	counter=\$((\$counter + 1))
	#done < ${bedfile}
	"""
}

process GroupReadsByUmi {
	tag "${Sample}"
	maxForks 5
	label 'process_medium'
	input:
		tuple val (Sample), path(mapped_bam)
	output: 
		tuple val (Sample), path ("*_grouped.bam"), emit : grouped_bam_ch
	script:
	"""
	for bam_file in ${mapped_bam}
	do 
		file_name=\$( basename \${bam_file} .bam)	# Removing the .bam extension
		java -Xmx${task.memory.toGiga()}g -jar ${params.fgbio_path} --compression 1 --async-io GroupReadsByUmi --input \${bam_file} --strategy Adjacency --edits 1 \
		--output \${file_name}_grouped.bam --family-size-histogram \${file_name}.tag-family-sizes.txt
	done
	# plot_family_sizes.py ${Sample}.tag-family-sizes.txt ${Sample}_family
	"""
}

process CallMolecularConsensusReads {
	tag "${Sample}"
	maxForks 5
	label 'process_medium'
	input:
		tuple val (Sample), path (grouped_bam)
	output:
		tuple val (Sample), path ("*_cons_umap.bam"), emit : consensus_bam_ch
	script:
	"""
	# This step generates unmapped consensus reads from the grouped reads and immediately filters them
	for bams in ${grouped_bam}
	do
		outfile_name=\$( basename \${bams} .bam)	# Removing the .bam extension
		java -Xmx${task.memory.toGiga()}g -jar ${params.fgbio_path} --compression 0 CallMolecularConsensusReads --input \${bams} \
		--output /dev/stdout --min-reads 4 --threads ${task.cpus} | java -Xmx${task.memory.toGiga()}g -jar ${params.fgbio_path} --compression 1 \
		FilterConsensusReads --input /dev/stdin --output \${outfile_name}_cons_umap.bam --ref ${params.genome} \
		--min-reads 4 --min-base-quality 20 --max-base-error-rate 0.25
	done
	"""
}

process COMBINEBAMS {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), path (cons_bam_list)
	output:
		tuple val (Sample), path ("${Sample}.bam"), emit : combined_bam_ch
	script:
	"""
	samtools merge -@ ${task.cpus} ${Sample}.bam ${cons_bam_list} 
	"""
}

process FilterConsBam {
	tag "${Sample}"
	label 'process_low'
	input:
		tuple val (Sample), file (cons_unmapped_bam)
	output:
		tuple val (Sample), file ("${Sample}_filt.bam"), file ("${Sample}_filt.bam.bai")
	script:
	"""
	${params.gatk} AddOrReplaceReadGroups I=${cons_unmapped_bam} O=${Sample}_with_rg.bam \
	RGID=AML RGLB=LIB-MIPS RGPL=ILLUMINA RGPU=UNIT_1 RGSM=${Sample}

	${params.samtools} fastq ${Sample}_with_rg.bam | bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200\\tPU:UNIT_1" -t ${task.cpus} -p -K 150000000 -Y ${params.genome} - | \
	java -Xmx${task.memory.toGiga()}g -jar ${params.fgbio_path} --compression 0 --async-io ZipperBams --unmapped ${cons_unmapped_bam} \
	--ref ${params.genome} --tags-to-reverse Consensus --tags-to-revcomp Consensus | ${params.samtools} sort --threads ${task.cpus} -o ${Sample}.cons.filtered.bam

	${params.gatk} AddOrReplaceReadGroups I=${Sample}.cons.filtered.bam O=${Sample}_filt.bam RGID=AML RGLB=LIB-MIPS RGPU=UNIT_1 RGPL=ILLUMINA RGSM=${Sample}
	${params.samtools} index ${Sample}_filt.bam > ${Sample}_filt.bam.bai
	"""
}

process SyntheticFastq {
	tag "${Sample}"
	input:
		tuple val (Sample), file (cons_filt_bam), file (cons_filt_bam_bai)
	output:
		tuple val (Sample), file ("*.synreads.bam"), file ("*.synreads.bam.bai")
	script:
	"""
	# Making a synthetic fastq and bam based on the number of reads in the sample bam file
	${params.PosControlScript} ${cons_filt_bam} ${Sample}_inreads

	#${params.fastq_bam} ${Sample}_inreads		# This script will output a file named *_inreads.fxd_sorted.bam

	bwa mem -R "@RG\\tID:AML\\tPL:ILLUMINA\\tLB:LIB-MIPS\\tSM:${Sample}\\tPI:200\\tPU:UNIT_1" -M -t ${task.cpus} ${params.genome} ${Sample}_inreads.fastq | \
	${params.samtools} sort -@ ${task.cpus} -o ${Sample}_inreads.fxd_sorted.bam -
	${params.samtools} index ${Sample}_inreads.fxd_sorted.bam > ${Sample}_inreads.fxd_sorted.bam.bai

	# Merging the Sample bam file with positive control bam file
	mv ${cons_filt_bam} ${Sample}.cons.filtered_merge.bam
	${params.samtools} merge -f ${Sample}_withsynreads.bam ${Sample}.cons.filtered_merge.bam ${Sample}_inreads.fxd_sorted.bam
	${params.gatk} AddOrReplaceReadGroups I=${Sample}_withsynreads.bam O=${Sample}.synreads.bam RGID=AML RGLB=LIB-MIPS RGPU=UNIT_1 RGPL=ILLUMINA RGSM=${Sample}
	${params.samtools} index ${Sample}.synreads.bam > ${Sample}.synreads.bam.bai
	"""
}