#!/usr/bin/bash
##BSUB -J smMIPS_pipeline
##BSUB -n 25
##BSUB -q normal
##-m cn2" to submit jobs to cn2
## or " -m cn3"

##########
#for ENTRY : BEDFILES#
##for LEUKEMIA/MIPS: /home/pipelines/mutation_detector_nextflow/bedfile/06112021_Leukemia_Panel_sorted
##for MIPS (IDT-MRD): /home/pipelines/mutation_detector_nextflow/bedfile/04243058_MRD_Panel_V1_final_sorted 
##for CNVpanel+ALP:/home/pipelines/mutation_detector_nextflow/bedfile/ALP_CNV_backbone_sorted
##for CNVpanel:/home/pipelines/mutation_detector_nextflow/bedfile/xgen-human-cnv-backbone-hyb-panel-probes
##for Lungpanel:/home/pipelines/mutation_detector_nextflow/bedfile/lung_panel_egfr_kras_tp53_sortd
##for Twistmyeloid:/home/pipelines/mutation_detector_nextflow/bedfile/Leukemia_Bed_file_MYFU_grch37_sorted
##for Twistlymphoid:/home/pipelines/mutation_detector_nextflow/bedfile/Leukemia_Bedfile_ALL_grch37hglft_genome_ucsc
##for combined_panel:/home/pipelines/mutation_detector_nextflow/bedfile/Leukemia_Panel_Myeloid_2023_Feb_hg37_sortd

echo "WARNING : change the bedfile and the cnv reference"
# for cnvkit reference 
# 06112021_Leukemia_Panel_sorted.bed : "/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_ref_GeneNames/Reference_labelled.cnn" 
# Leukemia_Panel_Myeloid_2023_Feb_hg37_sortd.bed : "/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/cnvkit_combpanel/Reference_combpanel.cnn"

source activate new_base
nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run main.nf -entry MRD \
--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd \
--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd \
--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
-resume -bg

#echo "This is outside the nextflow"

#for i in `cat samplesheet.csv`
#do 
#	cp Final_Output/$i/$i.final.cns Espresso_analysis/
#done

#cd Espresso_analysis
#./umi_error_model.sh /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv
#cd ../

conda deactivate 
