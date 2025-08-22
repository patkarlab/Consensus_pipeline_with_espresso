#!/bin/bash

#samplesheet=$1
#bedfile="/home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_060624_sortd.bed"

#for Sample in `cat ${samplesheet}`; do
#        /usr/bin/bedtools bamtobed -i Final_Output/${Sample}/${Sample}.uncollaps.bam > Final_Output/${Sample}/${Sample}.bed
#        /usr/bin/bedtools coverage -d -a ${bedfile} -b Final_Output/${Sample}/${Sample}.bed > Final_Output/${Sample}/${Sample}.uniform.bed
#        /home/pipelines/Consensus_pipeline_with_espresso/scripts/coverage_plots/coverage_plot.py Final_Output/${Sample}/${Sample}.uniform.bed Final_Output/${Sample}/${Sample}.uniform.bed.pdf
#done

#for Sample in `cat ${samplesheet}`; do
	#ls sequences/${Sample}*.fastq.gz 
	Sample=$1
	r1=sequences/${Sample}_S*_L*_R1_*.fastq.gz
	#r2=sequences/${Sample}_S*_L*_R2_*.fastq.gz
	base_count=$( zcat ${r1} | awk 'NR % 4 == 2 { bases += length($0);} END { print (bases*2)/10^9 }')
	echo "${Sample} ${base_count}"
#done		
