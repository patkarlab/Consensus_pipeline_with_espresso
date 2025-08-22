#!/usr/bin/bash

bamfile=$1
bedfile=$2
outdir=$3
Sample=$4

/usr/bin/bedtools bamtobed -i ${bamfile} > ${outdir}/temp.bed
#/usr/bin/bedtools coverage -d -a ${bedfile} -b ${outdir}/temp.bed > ${outdir}/${Sample}.uniform.bed

/usr/bin/bedtools coverage -counts -a ${bedfile} -b ${outdir}/temp.bed > ${outdir}/${Sample}_cov.regions_bedtools.bed

#/home/pipelines/Consensus_pipeline_with_espresso/scripts/coverage_plots/coverage_plot.py ${outdir}/${Sample}.uniform.bed ${outdir}/${Sample}_coverage_plot.pdf
