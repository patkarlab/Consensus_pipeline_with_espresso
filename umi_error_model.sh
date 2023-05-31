#!/usr/bin/bash

samplesheet=$1
bedfile=/home/pipelines/mutation_detector_nextflow/bedfile/04243058_MRD_Panel_V1_final_sorted.bed

#samtools sort /home/reference_genomes/NA12878/NA12878.final.bam > normal_sorted.bam
#samtools index normal_sorted.bam > normal_sorted.bam.bai
#bedtools sort -i ${bedfile} > sorted.bed

#source activate new_base
#for i in `cat ${samplesheet}`
#do 
	#echo ${i}
	#java -Xmx12g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 --async-io FastqToBam --input ${i}_*_R1_*.fastq.gz ${i}_*_R2_*.fastq.gz ${i}_*_R3_*.fastq.gz --read-structures +T +M +T --sample ${i} --library ${i} --output ${i}.unmapped.bam
	
	#samtools fastq ${i}.unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y /home/reference_genomes/hg19_broad/hg19_all.fasta - | java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 --async-io ZipperBams --unmapped ${i}.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --output ${i}.mapped.bam

	#java -Xmx8g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 --async-io GroupReadsByUmi --input ${i}.mapped.bam --strategy Adjacency --edits 1 --output ${i}.grouped.bam --family-size-histogram ${i}.tag-family-sizes.txt

	#java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 0 CallMolecularConsensusReads --input ${i}.grouped.bam --output /dev/stdout --min-reads 3 --min-input-base-quality 20 --threads 4 | java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 FilterConsensusReads --input /dev/stdin --output ${i}.cons.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --min-reads 3 --min-base-quality 45 --max-base-error-rate 0.2

#	samtools fastq ${i}.cons.unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y /home/reference_genomes/hg19_broad/hg19_all.fasta - | java -Xmx12g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 0 --async-io ZipperBams --unmapped ${i}.cons.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --tags-to-reverse Consensus --tags-to-revcomp Consensus | samtools sort --threads 8 -o ${i}.cons.filtered.bam --write-index
	
	# Making a synthetic fastq and bam based on the number of reads in the sample bam file 
#	./PosControlCount.sh ${i}.cons.filtered.bam ${i}_npm1
#	./fastq_bam.sh ${i}_npm1	# This script will output a file named *_npm1.fxd_sorted.bam

	# Merging the Sample bam file with positive control bam file
#	mv ${i}.cons.filtered.bam ${i}.cons.filtered_merge.bam
#	samtools merge -f ${i}.cons.filtered.bam ${i}.cons.filtered_merge.bam ${i}_npm1.fxd_sorted.bam
#	samtools index ${i}.cons.filtered.bam > ${i}.cons.filtered.bam.bai
	
	# ABRA2 for bamfile realignment
#	/usr/lib/jvm/java-8-openjdk-amd64/bin/java -Xmx16G -jar /home/programs/ABRA2/abra2-2.23.jar --in normal_sorted.bam,${i}.cons.filtered.bam --out normal.abra.bam,${i}.abra.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --threads 8 --targets sorted.bed --tmpdir ./ > abra.log

	# CNS file 
	#/home/programs/samtools-1.7/samtools mpileup -s -x -BQ 0 -q 1 -d 100000 --skip-indels -f /home/reference_genomes/hg19_broad/hg19_all.fasta ${i}.abra.bam > ${i}.mpileup
#	/home/programs/samtools-1.7/samtools mpileup -s -x -BQ 0 -q 1 -d 1000000 -f /home/reference_genomes/hg19_broad/hg19_all.fasta ${i}.abra.bam > ${i}.mpileup

	#/usr/lib/jvm/java-8-openjdk-amd64/bin/java -jar /home/programs/VarScan.v2.3.9.jar pileup2cns ${i}".mpileup" --variants SNP --min-coverage 10 --min-reads2 1 --min-avg-qual 30 --min-var-freq 0.0001 --p-value 1e-2 --strand-filter 0 > ${i}".cns"
#	/usr/lib/jvm/java-8-openjdk-amd64/bin/java -jar /home/programs/VarScan.v2.3.9.jar pileup2cns ${i}".mpileup" --min-coverage 10 --min-reads2 1 --min-avg-qual 30 --min-var-freq 0.0001 --p-value 1e-2 --strand-filter 0 > ${i}".cns"

#	grep -v '^chrM' ${i}".cns" > ${i}".nochrM.cns"
#	awk '{s=(NR==1)?"GT":"0/1";$0=$0 "\t" s}1' ${i}".nochrM.cns" > ${i}".final.cns"
#done
#source deactivate

#ESPRESSO
# R library
export R_LIBS="/home/tuhina/R/x86_64-pc-linux-gnu-library/3.6:$R_LIBS"
export R_LIBS="/usr/local/lib/R/site-library:$R_LIBS"
export R_LIBS="/usr/lib/R/site-library:$R_LIBS"
export R_LIBS="/usr/lib/R/library:$R_LIBS"

#cp /home/tuhina/MRD_DUPLEX/scripts/espresso_model.R $PWD/
#cp /home/tuhina/MRD_DUPLEX/nextflow_pipeline/variant_calling/COSMIC_heme_freq* $PWD/
#awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' ${bedfile} > mrd_bed_ranges.bed
#Rscript espresso_model.R

source activate new_base
perl /home/programs/annovar_latest/annovar/convert2annovar.pl -format vcf4 espresso_calls.vcf --outfile espresso.avinput --withzyg --includeinfo -allsample

for i in `cat ${samplesheet}`
do
	perl /home/programs/annovar_latest/annovar/table_annovar.pl "espresso.avinput."$i".avinput" --out $i --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all --operation g,r,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 /home/programs/annovar_latest/annovar/humandb/ --xreffile /home/programs/annovar_latest/annovar/example/gene_fullxref.txt

	python espresso_anno_format.py $i".hg19_multianno.csv" $i ${PWD}/

	cut -d',' -f1-28,47 $i"_espresso.csv" > $i"_final_anno.csv"
done

## MERGED EXCEL GENERATION
for i in `cat ${samplesheet}`
do
	python3 merge-csv_v3.py $i ./ $i.xlsx
done
