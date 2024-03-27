#!/usr/bin/bash

#bam_file=$1
mpileup_file=$1
min_base_qual=$2
outvcf=$3

source deactivate 
#samtools mpileup -F 0.2 -A -a -l /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd.bed -f "/home/reference_genomes/hg19_broad/hg19_all.fasta" -d 100000 ${bam_file} > temp.mpileup
perl /home/pipelines/Consensus_pipeline_with_espresso/scripts/call_variants_from_mpileup.pl ${mpileup_file} ${min_base_qual} ${outvcf}

#perl /home/programs/annovar_latest/annovar/convert2annovar.pl -format vcf4 ${outvcf}  --outfile temp.avinput --withzyg --includeinfo

#perl /home/programs/annovar_latest/annovar/table_annovar.pl temp.avinput --out 24NGS227-MRD.mips --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 /home/programs/annovar_latest/annovar/humandb/ --xreffile /home/programs/annovar_latest/annovar/example/gene_fullxref.txt


