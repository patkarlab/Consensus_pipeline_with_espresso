#!/bin/bash

samplesheet=$1
bedfile="/home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_060624_sortd.bed"

for Sample in `cat ${samplesheet}`; do
        /usr/bin/bedtools bamtobed -i Final_Output/${Sample}/${Sample}.uncollaps.bam > Final_Output/${Sample}/${Sample}.bed
        /usr/bin/bedtools coverage -d -a ${bedfile} -b Final_Output/${Sample}/${Sample}.bed > Final_Output/${Sample}/${Sample}.uniform.bed
        /home/pipelines/Consensus_pipeline_with_espresso/scripts/coverage_plots/coverage_plot.py Final_Output/${Sample}/${Sample}.uniform.bed Final_Output/${Sample}/${Sample}.uniform.bed.pdf
done

