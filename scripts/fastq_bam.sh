#!/usr/bin/bash

input_basename=$1

bwa mem -R "@RG\tID:AML\tPL:ILLUMINA\tLB:LIB-MIPS\tSM:${input_basename}\tPI:200" -M -t 20 /home/reference_genomes/hg19_broad/hg19_all.fasta ${input_basename}.fastq > ${input_basename}.sam

/home/programs/samtools-1.7/samtools view -bT /home/reference_genomes/hg19_broad/hg19_all.fasta ${input_basename}.sam > ${input_basename}.fxd.bam
/home/programs/samtools-1.7/samtools sort ${input_basename}.fxd.bam > ${input_basename}.fxd_sorted.bam
/home/programs/samtools-1.7/samtools index ${input_basename}.fxd_sorted.bam > ${input_basename}.fxd_sorted.bam.bai
