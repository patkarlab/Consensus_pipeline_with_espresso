#!/usr/bin/bash

#samplesheet=$1
mutect2=$1
varscan=$2
vardict=$3

#for i in `cat ${samplesheet}` 
#do 
perl /home/programs/annovar_latest/annovar//convert2annovar.pl -format vcf4 ${mutect2} --outfile mutect2.avinput --withzyg --includeinfo	
perl /home/programs/annovar_latest/annovar//table_annovar.pl mutect2.avinput --out mutect2.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 /home/programs/annovar_latest/annovar//humandb/ --xreffile /home/programs/annovar_latest/annovar//example/gene_fullxref.txt
	
if [ -s mutect2.somaticseq.hg19_multianno.csv ]; then
	#python3 /home/pipelines/Consensus_pipeline_with_espresso/scripts/somaticseqoutput-format_v2.py NPM1-078.somaticseq.hg19_multianno.csv ${i}.csv
	ls mutect2.somaticseq.hg19_multianno.csv
else
	touch mutect2.somaticseq.hg19_multianno.csv
fi
#done

perl /home/programs/annovar_latest/annovar//convert2annovar.pl -format vcf4 ${varscan} --outfile varscan.avinput --withzyg --includeinfo
perl /home/programs/annovar_latest/annovar//table_annovar.pl varscan.avinput --out varscan.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 /home/programs/annovar_latest/annovar//humandb/ --xreffile /home/programs/annovar_latest/annovar//example/gene_fullxref.txt

if [ -s varscan.somaticseq.hg19_multianno.csv ]; then
	ls varscan.somaticseq.hg19_multianno.csv
else 
	touch varscan.somaticseq.hg19_multianno.csv
fi


perl /home/programs/annovar_latest/annovar//convert2annovar.pl -format vcf4 ${vardict} --outfile vardict.avinput --withzyg --includeinfo
perl /home/programs/annovar_latest/annovar//table_annovar.pl vardict.avinput --out vardict.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 /home/programs/annovar_latest/annovar//humandb/ --xreffile /home/programs/annovar_latest/annovar//example/gene_fullxref.txt

if [ -s vardict.somaticseq.hg19_multianno.csv ]; then
	ls vardict.somaticseq.hg19_multianno.csv
else 
	touch vardict.somaticseq.hg19_multianno.csv
fi	
