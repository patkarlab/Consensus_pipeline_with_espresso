#!/bin/bash

Sample=$1

#Trimming
trimmomatic PE  /home/pipelines/Consensus_pipeline_with_espresso/sequences//${Sample}_*R1_*.fastq.gz /home/pipelines/spresso/sequences//${Sample}_*R2_*.fastq.gz       -baseout ${Sample}.fq.gz  ILLUMINACLIP:/home/programs/Trimmomatic-0.32:30:10:2:keepBothReads         ILLUMINACLIP:/home/programs/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:keeIDINGWINDOW:4:15 MINLEN:40
sleep 5s

#FastqToBam
java -Xmx12g -jar "/home/programs/fgbio/fgbio-2.0.1.jar" --compression 1 --async-io FastqToBam --input ${Sample}_1P.fq.gz ${Sample}_2P.fq.gz --read-structures 4M+T 4M+T --sample ${Sample} --library ${Sample} --output ${Sample}.unmapped.bam

#MapBam
/home/programs/samtools-1.7/samtools fastq ${Sample}.unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y /home/reference_genomes/hg19_broad/hg19_all.fasta - | java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 --async-io ZipperBams --unmapped ${Sample}.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --output ${Sample}.mapped.bam

#GroupReadsByUmi
java -Xmx8g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 --async-io GroupReadsByUmi --input ${Sample}.mapped.bam --strategy Adjacency --edits 1 --output ${Sample}.grouped.bam --family-size-histogram ${Sample}.tag-family-sizes.txt

#CallMolecularConsensusReads
java -Xmx25g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 0 CallMolecularConsensusReads --input ${Sample}.grouped.bam --output /dev/stdout --min-reads 4 --threads 4 | java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 FilterConsensusReads --input /dev/stdin --output ${Sample}.cons.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --min-reads 4 --min-base-quality 20 --max-base-error-rate 0.25

#FilterConsBam
/home/programs/samtools-1.7/samtools fastq ${Sample}.cons.unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y /home/reference_genomes/hg19_broad/hg19_all.fasta - | java -Xmx12g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 0 --async-io ZipperBams --unmapped ${Sample}.cons.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --tags-to-reverse Consensus --tags-to-revcomp Consensus | /home/programs/samtools-1.7/samtools sort --threads 8 -o ${Sample}.cons.filtered.bam


/home/programs/samtools-1.7/samtools index ${Sample}.cons.filtered.bam > ${Sample}.cons.filtered.bam.bai


#ABRA2_realign
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -Xmx16G -jar /home/programs/ABRA2/abra2-2.23.jar --in ${Sample}.cons.filtered.bam --out ${Sample}.ErrorCorrectd.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --threads 8 --targets /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_asxl2ex11_sortd.bed --tmpdir ./ > abra.log


/home/programs/samtools-1.7/samtools index ${Sample}.ErrorCorrectd.bam > ${Sample}.ErrorCorrectd.bam.bai


#Coverage
/home/pipelines/Consensus_pipeline_with_espresso/scripts/mosdepth.sh ${Sample}.ErrorCorrectd.bam ${Sample}_umi_cov /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_asxl2ex11_sortd.bed


#hsmetrics
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -jar /home/programs/picard/build/libs/picard.jar CollectHsMetrics I= ${Sample}.ErrorCorrectd.bam O= ${Sample}_hsmetrics.txt BAIT_INTERVALS= /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_asxl2ex11_sortd.interval_list TARGET_INTERVALS= /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_asxl2ex11_sortd.interval_list R= /home/reference_genomes/hg19_broad/hg19_all.fasta VALIDATION_STRINGENCY=LENIENT


#Mutect2
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -Xmx10G -jar /home/diagnostics/programs/gatk-4.2.6.0/gatk-package-4.2.6.0-local.jar Mutect2 -R /home/reference_genomes/hg19_broad/hg19_all.fasta -I:tumor ${Sample}.ErrorCorrectd.bam -O ${Sample}.mutect2.vcf -L /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_asxl2ex11_sortd.bed -mbq 20


/usr/lib/jvm/java-8-openjdk-amd64/bin//java -Xmx10G -jar /home/diagnostics/programs/gatk-4.2.6.0/gatk-package-4.2.6.0-local.jar FilterMutectCalls -V ${Sample}.mutect2.vcf --stats ${Sample}.mutect2.vcf.stats -O ${Sample}_filtered.vcf -R /home/reference_genomes/hg19_broad/hg19_all.fasta --unique-alt-read-count 3 --min-median-base-quality 20 --min-median-mapping-quality 30

#VarDict
VarDict -G /home/reference_genomes/hg19_broad/hg19_all.fasta -f 0.0001 -r 8 -N ${Sample} -b ${Sample}.ErrorCorrectd.bam -c 1 -S 2 -E 3 -g 4 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/mips_mrd_bal11March_asxl2ex11_sortd.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N ${Sample} -E -f 0.0001 > ${Sample}.vardict.vcf


#Varscan
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -jar /home/programs/VarScan.v2.3.9.jar mpileup2indel ${Sample}.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 1 --min-var-freq 0.000001 --p-value 1e-1 --output-vcf 1 > ${Sample}.varscan_indel.vcf
bgzip -c ${Sample}.varscan_snp.vcf > ${Sample}.varscan_snp.vcf.gz
bgzip -c ${Sample}.varscan_indel.vcf > ${Sample}.varscan_indel.vcf.gz
/home/programs/bcftools-1.9/bcftools index -t ${Sample}.varscan_snp.vcf.gz
/home/programs/bcftools-1.9/bcftools index -t ${Sample}.varscan_indel.vcf.gz
/home/programs/bcftools-1.9/bcftools concat -a ${Sample}.varscan_snp.vcf.gz ${Sample}.varscan_indel.vcf.gz -o ${Sample}.varscan.vcf
