# fgbio concensus pipeline with Espresso

## Description
This pipeline follows the fgbio pipeline to obtain mapped consensus bam file. This bam file is used as input to VarScan, the output of which is used to perform error corrected variant calling using Espresso

## 1. Demultiplexing the umi tagged sequences
### Use the following command if index1 has a 8 base umi at the end
```
bcl2fastq -R 230215_M04898_0087_000000000-DJG93 -o 230215_M04898_0087_000000000-DJG93/demux_out/ --sample-sheet Samplesheet_delN_16Feb2023.csv --use-bases-mask Y145,I8,Y8,I8,Y153 --mask-short-adapter-reads 0
```
Here, -R is the input (entire run) folder , -o is the output folder and --sample-sheet is the miseq format sample sheet (i2 is reverse complimented). The original samplesheet used was SampleSheetUsed.csv.  

## 2. Obtaining consensus bam file
The umi_error_model.sh follows the fgbio pipeline mentioned here https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md to obtain a consensus bam file. This script requires a normal bam file, a bed file, samplesheet.csv.
It can be run using 
```
./umi_error_model.sh samplesheet.csv
```
This script uses PosControlCount.sh which adds synthetic npm1 reads.  
It requires MRD_Positive_Control_wt.fastq and MRD_Positive_Control.fastq which have wild type and insert containing sequence respectively.  
This script counts the total no. of reads in the input bam file. The number of npm1 reads to be added in the fastq file are calculated such that they form 6% of the total reads in the final merged bam file. The npm1 insert has a final VAF of 4%.  
The fastq_bam.sh script converts the fastq file obtained from the PosControlCount.sh to bam file. This bam file is then integrated in the sample's bam file.  

## 3. Varscan for variant calling
The bam file from the above step is used to obtain a mpileup file. 
This mpileup file is used as input for VarScan for variant calling.  

## 4. Espresso for error corrected variant detection
The .cns files obtained from VarScan are provided as input to Espresso.  The output from Espresso is a single vcf file for all the samples.  
The detailed espresso workflow is present here : https://htmlpreview.github.io/?https://github.com/abelson-lab/Espresso/blob/master/vignettes/Espresso_workflow.html  
The github repo for espresso : https://github.com/abelson-lab/Espresso  

## 5. Annotation of variants
The .vcf file is separated as per the samples and annotated using ANNOVAR.  
The annotated csv files are converted to excel files.  
