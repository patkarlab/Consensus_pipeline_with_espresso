#!/usr/bin/bash
##BSUB -J smMIPS_pipeline
##BSUB -n 25
##BSUB -q normal
##-m cn2" to submit jobs to cn2
## or " -m cn3"

##########
# Original bedfiles 
# bedfile: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd, bedfile2: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd

# Updated bedfiles 22Dec23
# bedfile: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_updated_uniq_sortd, bedfile2: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_updated_sortd

source activate new_base
nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run main.nf -entry MRD \
--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd \
--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd \
--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
-resume -bg
conda deactivate 
