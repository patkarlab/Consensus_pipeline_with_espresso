#!/usr/bin/bash
##########
# Original bedfiles 
# bedfile: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd, bedfile2: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd

# Updated bedfiles 22Dec23
# bedfile: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_updated_uniq_sortd, bedfile2: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_updated_sortd

# Updated bedfile 24Jan24
# bedfile: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_24jan24_sortd, bedfile2: /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd

# Original bedfile for CEBPA-MRD
#bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_9857C1725CE84D37803C380036834C74_g_sortd, bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd

# For AML-MRD
#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run main.nf -entry MRD \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_24jan24_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd \
#--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
#--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
#-resume -bg
#conda deactivate 

#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run mips_mrd.nf -entry MRD \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal15March_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal15March_exon_sortd \
#--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
#--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
#-resume -bg
#conda deactivate

# For CEBPA-MRD
#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run main.nf -entry MRD \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/10May24_CEBPA_MRD/IDT_MRD_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd \
#--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
#--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
#-resume -bg
#conda deactivate

#For MIPS-MRD
#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run mips_mrd.nf -entry MRD \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_asxl2ex11_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_exon_sortd \
#--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
#--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
#-resume -bg
#conda deactivate

# For Narasimha MRD
source activate new_base
nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run mips_mrd.nf -entry NARASIMHA_MRD \
--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_060624_hotspot_sortd \
--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_exon_sortd \
--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
-resume -bg
conda deactivate

# For Narasimha MRD Dragen
#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run narasimha_mrd_dragen.nf -entry NARASIMHA_MRD \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_060624_hotspot_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_exon_sortd \
#--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
#--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
#-resume -bg
#conda deactivate
