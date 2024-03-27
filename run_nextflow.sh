#!/usr/bin/bash
##########
# Original bedfiles 
# bedfile: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd, bedfile2: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd

# Updated bedfiles 22Dec23
# bedfile: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_updated_uniq_sortd, bedfile2: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_updated_sortd

# Updated bedfile 24Jan24
# bedfile: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_24jan24_sortd, bedfile2: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd

#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run main.nf -entry MRD \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_24jan24_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd \
#--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
#--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
#-resume -bg
#conda deactivate 

source activate new_base
nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run mips_mrd.nf -entry MRD \
--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal15March_sortd \
--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal15March_exon_sortd \
--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
-resume -bg
conda deactivate
