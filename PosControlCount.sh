#! /usr/bin/env bash
#This script will take a bam file as input, count the no. of reads and make a fastq file with a 4% npm1 VAF

infile=$1
outfile=$2

npm1_reads=$(samtools view -c ${infile} | awk '{print int ($1*(6/94))}')  # No. of extra reads to be added, 6% of the final reads
touch ${outfile}.fastq
ref=$(echo ${npm1_reads} | awk '{print int ($1*0.96)}')	# 96% of the reads are ref
alt=$(echo ${npm1_reads} | awk '{print int ($1*0.04)}')	# 4% of the reads are alt
for ((j=1; j<=${ref}; j++)); do cat MRD_Positive_Control_wt.fastq; done > ${outfile}.fastq
for ((k=1; k<=${alt}; k++)); do cat MRD_Positive_Control.fastq; done >> ${outfile}.fastq
