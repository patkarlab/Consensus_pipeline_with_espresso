# fgbio concensus pipeline with Espresso

## Description
This pipeline follows the fgbio pipeline to obtain mapped consensus bam file. This bam file is used as input to VarScan, the output of which is used to perform error corrected variant calling using Espresso

## 1. Demultiplexing the umi tagged sequences
## Use the following command if index1 has a 8 base umi at the end
```
bcl2fastq -R 230215_M04898_0087_000000000-DJG93 -o 230215_M04898_0087_000000000-DJG93/demux_out/ --sample-sheet Samplesheet_delN_16Feb2023.csv --use-bases-mask Y145,I8,Y8,I8,Y153 --mask-short-adapter-reads 0
```

## 2. Obtaining consensus bam file
The umi_error_model.sh follows the fgbio pipeline mentioned here https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md to obtain a consensus bam file. This script requires a normal bam file, a bed file, samplesheet.csv
It can be run using 
```
./umi_error_model.sh samplesheet.csv
```
