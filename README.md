# MIPS pipeline

## Description
This is a nextflow pipeline for analysing data for MIPS-MRD. It follows the fgbio pipeline (https://github.com/fulcrumgenomics/fgbio) to obtain mapped consensus bam file. 

For running this pipeline, following programs need to be installed and their complete paths need to be added in the params section of the `nextflow.config`.

- fastqc = fastqc executable path 
- java_path = directory containing the java executable
- GATK38_path = path to the GenomeAnalysisTK-3.8.jar file
- GATK42_path = path to the gatk-package-4.2.6.0-local.jar file
- fgbio_path = ath to the fgbio-2.0.1.jar file
- samtools = samtools executable path
- genome = Genomic fasta file
- site1 = known_polymorphic_sites 1 (Mills_and_1000G_gold_standard.indels)
- PosControlScript = Custom script for generating Synthetic Fastq with 4% npm1 VAF
- fastq_bam = Custom script for generating bam from fastq with bwa
- abra2_path = directory containing the abra jar file
- pear_path = pear executable path 
- bedtools = bedtools executable path
- mosdepth
- picard_path = path to the picard.jar file
- get_itd_path = path to the folder containing getitd.py
- VarDict
- varscan_path = path to the VarScan.v2.3.9.jar file
- annovarLatest_path = path to the folder containing ANNOVAT perl scripts
- bcftools_path = bcftools executable path 
- somaticseq

## Usage:

1. Keep the `fastq` files into the `sequences/` folder.

2. Change the `samplesheet.csv`. It should have a list of IDs of the samples. 

3.  Run the following script.

```
./run_nextflow.sh > script.log
```
This script contains the nextflow command used to execute the workflow.

```
source activate new_base

nextflow -c /home/pipelines/Consensus_pipeline_with_espresso/nextflow.config run mips_mrd.nf -entry MRD \
--bedfile /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal210125_sortd \
--bedfile2 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_exon_sortd \
--sequences /home/pipelines/Consensus_pipeline_with_espresso/sequences/ \
--input /home/pipelines/Consensus_pipeline_with_espresso/samplesheet.csv \
-resume -bg

conda deactivate
```
