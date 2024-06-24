#!/usr/bin/bash

/usr/bin/bedtools bamtobed -i 24NGS631-NV2-AMLMRD.uncollaps.bam > 24NGS631-NV2-AMLMRD.bed
/usr/bin/bedtools coverage -d -a /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_060624_sortd.bed -b 24NGS631-NV2-AMLMRD.bed > 24NGS631-NV2-AMLMRD.uniform.bed
./coverage_plot.py 24NGS631-NV2-AMLMRD.uniform.bed temp.pdf
