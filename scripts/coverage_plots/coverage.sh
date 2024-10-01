#!/usr/bin/bash

bamfile=$1
bedfile=$2
outdir=$3

/usr/bin/bedtools bamtobed -i ${bamfile} > ${outdir}/temp.bed
/usr/bin/bedtools coverage -d -a ${bedfile} -b ${outdir}/temp.bed > ${outdir}/temp.uniform.bed
/home/pipelines/Consensus_pipeline_with_espresso/scripts/coverage_plots/coverage_plot.py ${outdir}/temp.uniform.bed ${outdir}/coverage_plot.pdf
rm ${outdir}/temp.bed ${outdir}/temp.uniform.bed
