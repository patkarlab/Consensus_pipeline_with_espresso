#!/usr/bin/bash

samplesheet=$1
bedfile=$2
#/home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed

#ESPRESSO
# R library
export R_LIBS="/home/tuhina/R/x86_64-pc-linux-gnu-library/3.6:$R_LIBS"
#export R_LIBS="/usr/local/lib/R/site-library:$R_LIBS"
#export R_LIBS="/usr/lib/R/site-library:$R_LIBS"
#export R_LIBS="/usr/lib/R/library:$R_LIBS"


awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' ${bedfile} > mrd_bed_ranges.bed
#Rscript espresso_model.R
source deactivate

Rscript /home/pipelines/Consensus_pipeline_with_espresso/scripts/Espresso/espresso_model.R

source activate new_base
perl /home/programs/annovar_latest/annovar/convert2annovar.pl -format vcf4 espresso_calls.vcf --outfile espresso.avinput --withzyg --includeinfo -allsample

for i in `cat ${samplesheet}`
do
	perl /home/programs/annovar_latest/annovar/table_annovar.pl "espresso.avinput."$i".final.cns.avinput" --out $i --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all --operation g,r,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 /home/programs/annovar_latest/annovar/humandb/ --xreffile /home/programs/annovar_latest/annovar/example/gene_fullxref.txt

	python /home/pipelines/Consensus_pipeline_with_espresso/scripts/espresso_anno_format.py $i".hg19_multianno.csv" $i ${PWD}/
	cut -d',' -f1-28,47 $i"_espresso.csv" > $i"_final_anno.csv"
done

## MERGED EXCEL GENERATION
for i in `cat ${samplesheet}`
do
	python3 /home/pipelines/Consensus_pipeline_with_espresso/scripts/merge-csv_v3.py $i ./ /home/pipelines/Consensus_pipeline_with_espresso/Final_Output/${i}/${i}.Espresso.xlsx
done
