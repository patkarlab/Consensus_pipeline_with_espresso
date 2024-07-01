#!/bin/bash

for Sample in `cat ../samplesheet.csv`; do
    /usr/bin/bedtools bamtobed -i "${Sample}/${Sample}.uncollaps.bam" > "${Sample}/${Sample}.uncollaps.bed"
    /usr/bin/bedtools coverage -counts -a /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_060624_sortd.bed -b "${Sample}/${Sample}.uncollaps.bed" > "${Sample}/${Sample}_uncollaps.regions_bedtools.bed"

	/usr/bin/bedtools bamtobed -i ${Sample}/${Sample}.ErrorCorrectd.bam > ${Sample}/${Sample}.ErrorCorrectd.bed
	/usr/bin/bedtools coverage -counts -a /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_060624_sortd.bed -b ${Sample}/${Sample}.ErrorCorrectd.bed > ${Sample}/${Sample}.ErrorCorrectd_bedtools.bed
done
