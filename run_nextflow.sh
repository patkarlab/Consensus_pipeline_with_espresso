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
source activate new_base
nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run scripts/main.nf -entry MRD \
--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_IDT_192_withoutCEBPA_sortd \
--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_withoutCEBPA_sortd \
--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
--input /home/pipelines/Consensus_pipeline_with_espresso/temp_samp.dat \
--outdir /home/pipelines/Consensus_pipeline_with_espresso/Final_Output \
-resume -bg
conda deactivate 

# For old MIPS-MRD
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
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run scripts/main.nf -entry MRD \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/10May24_CEBPA_MRD/IDT_MRD_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd \
#--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
#--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
#-resume -bg
#conda deactivate

#For MIPS-MRD
#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run mips_mrd.nf -entry MRD \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal210125_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_exon_sortd \
#--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
#--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
#--outdir /home/pipelines/Consensus_pipeline_with_espresso/Final_Output \
#-resume -bg
#conda deactivate

# For Narasimha MRD
#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run mips_mrd.nf -entry NARASIMHA_MRD \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal210125_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_exon_sortd \
#--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
#--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
#-resume -bg
#conda deactivate

# For Narasimha MRD Dragen
#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run narasimha_mrd_dragen.nf -entry NARASIMHA_MRD \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/CEBPA_MRD_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/CEBPA_MRD_sortd \
#--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
#--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
#-resume -bg
#conda deactivate

# MIPS-gencore trial
#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run scripts/mips_gencore.nf -entry MIPS_MRD_Gencore \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_asxl2ex11_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_exon_sortd \
#-resume -bg

# NPM1-FLT3 MRD
#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/scripts/npm1_flt3.config run scripts/npm1_flt3_mrd.nf -entry NPM1 \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/NPM1_FLT3 \
#-resume -bg

# For Narasimha MRD
# source activate new_base
# nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run mips_mrd.nf -entry U2AF1MRD \
# --bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/u2af1_sortd \
# --bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/u2af1_sortd \
# -resume -bg 
# conda deactivate

# For U2AF1 MRD
#source activate new_base
#nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run U2AF1-MRD.nf -entry U2AF1MRD_FGBIO \
#--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/CEBPA_MRD_sortd \
#--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/CEBPA_MRD_sortd \
#-resume -bg 
#conda deactivate
