#!/usr/bin/bash

source activate new_base

# FastqToBam
java -Xmx12g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 --async-io FastqToBam --input /home/pipelines/Consensus_pipeline_with_espresso/sequences//NPM1-078_*_R1_*.fastq.gz /home/pipelines/Consensus_pipeline_with_espresso/sequences//NPM1-078_*_R2_*.fastq.gz /home/pipelines/Consensus_pipeline_with_espresso/sequences//NPM1-078_*_R3_*.fastq.gz --read-structures +T +M +T --sample NPM1-078 --library NPM1-078 --output NPM1-078.unmapped.bam

# MapBam
/home/programs/samtools-1.7/samtools fastq NPM1-078.unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y /home/reference_genomes/hg19_broad/hg19_all.fasta - | java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 --async-io ZipperBams --unmapped NPM1-078.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --output NPM1-078.mapped.bam

/home/programs/samtools-1.7/samtools sort NPM1-078.mapped.bam -o NPM1-078.Uncollaps.bam
/home/programs/samtools-1.7/samtools index NPM1-078.Uncollaps.bam > NPM1-078.Uncollaps.bam.bai

# GroupReadsByUmi
java -Xmx8g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 --async-io GroupReadsByUmi --input NPM1-078.Uncollaps.bam --strategy Adjacency --edits 1 --output NPM1-078.grouped.bam --family-size-histogram NPM1-078.tag-family-sizes.txt

# vardict_uncoll
VarDict -G /home/reference_genomes/hg19_broad/hg19_all.fasta -f 0.0001 -N NPM1-078 -b NPM1-078.Uncollaps.bam -c 1 -S 2 -E 3 -g 4 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N NPM1-078 -E -f 0.0001 > NPM1-078.vardict.vcf

# hsmetrics_run_uncoll
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -jar /home/programs/picard/build/libs/picard.jar CollectHsMetrics I= NPM1-078.Uncollaps.bam O= NPM1-078_uncoll_hsmetrics.txt BAIT_INTERVALS= /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.interval_list TARGET_INTERVALS= /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.interval_list R= /home/reference_genomes/hg19_broad/hg19_all.fasta VALIDATION_STRINGENCY=LENIENT

# varscan_uncoll
/home/programs/samtools-1.7/samtools mpileup -d 1000000 -f /home/reference_genomes/hg19_broad/hg19_all.fasta NPM1-078.Uncollaps.bam > NPM1-078.mpileup
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -jar /home/programs/VarScan.v2.3.9.jar mpileup2snp NPM1-078.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1e-4 --output-vcf 1 > NPM1-078.varscan_snp.vcf
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -jar /home/programs/VarScan.v2.3.9.jar mpileup2indel NPM1-078.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1e-4 --output-vcf 1 > NPM1-078.varscan_indel.vcf
bgzip -c NPM1-078.varscan_snp.vcf > NPM1-078.varscan_snp.vcf.gz
bgzip -c NPM1-078.varscan_indel.vcf > NPM1-078.varscan_indel.vcf.gz
/home/programs/bcftools-1.9/bcftools index -t NPM1-078.varscan_snp.vcf.gz
/home/programs/bcftools-1.9/bcftools index -t NPM1-078.varscan_indel.vcf.gz
/home/programs/bcftools-1.9/bcftools concat -a NPM1-078.varscan_snp.vcf.gz NPM1-078.varscan_indel.vcf.gz -o NPM1-078.varscan.vcf

# mutect2_run_uncoll
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -Xmx10G -jar /home/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T MuTect2 -R /home/reference_genomes/hg19_broad/hg19_all.fasta -I:tumor NPM1-078.Uncollaps.bam -o NPM1-078.mutect2.vcf --dbsnp /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.vcf -L /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed -nct 25 -contamination 0.02 -mbq 20

# coverage_mosdepth_uncoll
/home/pipelines/Consensus_pipeline_with_espresso/scripts/mosdepth.sh NPM1-078.Uncollaps.bam NPM1-078_uncollaps /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed
/home/pipelines/Consensus_pipeline_with_espresso/scripts/mosdepth.sh NPM1-078.Uncollaps.bam NPM1-078_exon_uncollaps /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd.bed

# lofreq_uncoll
/home/programs/lofreq_star-2.1.4_linux-x86-64/bin/lofreq viterbi -f /home/reference_genomes/hg19_broad/hg19_all.fasta -o NPM1-078.lofreq.pre.bam NPM1-078.Uncollaps.bam
/home/programs/samtools-1.7/samtools sort NPM1-078.lofreq.pre.bam > NPM1-078.lofreq.bam
/home/programs/lofreq_star-2.1.4_linux-x86-64/bin/lofreq call -b dynamic -C 2 -a 0.00005 -q 20 -Q 20 -m 50 -f /home/reference_genomes/hg19_broad/hg19_all.fasta -l /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed -o NPM1-078.lofreq.vcf NPM1-078.lofreq.bam --no-default-filter
/home/programs/lofreq_star-2.1.4_linux-x86-64/bin/lofreq filter --no-defaults -a 0.0001 -i NPM1-078.lofreq.vcf -o NPM1-078.lofreq.filtered.vcf

# somaticSeq_run_uncoll
somaticseq_parallel.py --output-directory NPM1-078.somaticseq --genome-reference /home/reference_genomes/hg19_broad/hg19_all.fasta --inclusion-region /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed --threads 25 --algorithm xgboost  --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file NPM1-078.Uncollaps.bam --mutect2-vcf NPM1-078.mutect2.vcf --vardict-vcf NPM1-078.vardict.vcf --varscan-vcf NPM1-078.varscan.vcf --lofreq-vcf NPM1-078.lofreq.filtered.vcf --sample-name NPM1-078

/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/vcf_sorter.sh NPM1-078.somaticseq/Consensus.sSNV.vcf NPM1-078.somaticseq/somaticseq_snv.vcf
bgzip -c NPM1-078.somaticseq/somaticseq_snv.vcf > NPM1-078.somaticseq/somaticseq_snv.vcf.gz
/home/programs/bcftools-1.9/bcftools index -t NPM1-078.somaticseq/somaticseq_snv.vcf.gz

/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/vcf_sorter.sh NPM1-078.somaticseq/Consensus.sINDEL.vcf NPM1-078.somaticseq/somaticseq_indel.vcf
bgzip -c NPM1-078.somaticseq/somaticseq_indel.vcf > NPM1-078.somaticseq/somaticseq_indel.vcf.gz
/home/programs/bcftools-1.9/bcftools index -t NPM1-078.somaticseq/somaticseq_indel.vcf.gz

/home/programs/bcftools-1.9/bcftools concat -a NPM1-078.somaticseq/somaticseq_snv.vcf.gz NPM1-078.somaticseq/somaticseq_indel.vcf.gz -o NPM1-078.UncollapsVariant.vcf

perl /home/programs/annovar_latest/annovar//convert2annovar.pl -format vcf4 NPM1-078.UncollapsVariant.vcf --outfile NPM1-078.UncollapsVariant.avinput --withzyg --includeinfo

perl /home/programs/annovar_latest/annovar//table_annovar.pl NPM1-078.UncollapsVariant.avinput --out NPM1-078.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 /home/programs/annovar_latest/annovar//humandb/ --xreffile /home/programs/annovar_latest/annovar//example/gene_fullxref.txt

if [ -s NPM1-078.somaticseq.hg19_multianno.csv ]; then
        python3 /home/pipelines/Consensus_pipeline_with_espresso/scripts/somaticseqoutput-format_v2.py NPM1-078.somaticseq.hg19_multianno.csv NPM1-078.UncollapsVariant.csv
else
        touch NPM1-078.UncollapsVariant.csv
fi

# CallMolecularConsensusReads
java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 0 CallMolecularConsensusReads --input NPM1-078.grouped.bam --output /dev/stdout --min-reads 3 --min-input-base-quality 20 --threads 4 | java -Xmx4g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 1 FilterConsensusReads --input /dev/stdin --output NPM1-078.cons.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --min-reads 3 --min-base-quality 45 --max-base-error-rate 0.2

# FilterConsBam
/home/programs/samtools-1.7/samtools fastq NPM1-078.cons.unmapped.bam | bwa mem -t 16 -p -K 150000000 -Y /home/reference_genomes/hg19_broad/hg19_all.fasta - | java -Xmx12g -jar /home/programs/fgbio/fgbio-2.0.1.jar --compression 0 --async-io ZipperBams --unmapped NPM1-078.cons.unmapped.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --tags-to-reverse Consensus --tags-to-revcomp Consensus | /home/programs/samtools-1.7/samtools sort --threads 8 -o NPM1-078.cons.filtered.bam

# SyntheticFastq
# Making a synthetic fastq and bam based on the number of reads in the sample bam file
/home/pipelines/Consensus_pipeline_with_espresso/scripts/PosControlCount.sh NPM1-078.cons.filtered.bam NPM1-078_inreads
/home/pipelines/Consensus_pipeline_with_espresso/scripts/fastq_bam.sh NPM1-078_inreads          # This script will output a file named *_inreads.fxd_sorted.bam

# Merging the Sample bam file with positive control bam file
mv NPM1-078.cons.filtered.bam NPM1-078.cons.filtered_merge.bam
/home/programs/samtools-1.7/samtools merge -f NPM1-078.synreads.bam NPM1-078.cons.filtered_merge.bam NPM1-078_inreads.fxd_sorted.bam
/home/programs/samtools-1.7/samtools index NPM1-078.synreads.bam > NPM1-078.synreads.bam.bai

# ABRA2_realign
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -Xmx16G -jar /home/programs/ABRA2/abra2-2.23.jar --in NPM1-078.synreads.bam --out NPM1-078.ErrorCorrectd.bam --ref /home/reference_genomes/hg19_broad/hg19_all.fasta --threads 8 --targets /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed --tmpdir ./ > abra.log
/home/programs/samtools-1.7/samtools index NPM1-078.ErrorCorrectd.bam > NPM1-078.ErrorCorrectd.bam.bai

# vardict
VarDict -G /home/reference_genomes/hg19_broad/hg19_all.fasta -f 0.0001 -N NPM1-078 -b NPM1-078.ErrorCorrectd.bam -c 1 -S 2 -E 3 -g 4 /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed | sed '1d' | teststrandbias.R | var2vcf_valid.pl -N NPM1-078 -E -f 0.0001 > NPM1-078.vardict.vcf

# lofreq
/home/programs/lofreq_star-2.1.4_linux-x86-64/bin/lofreq viterbi -f /home/reference_genomes/hg19_broad/hg19_all.fasta -o NPM1-078.lofreq.pre.bam NPM1-078.ErrorCorrectd.bam
/home/programs/samtools-1.7/samtools sort NPM1-078.lofreq.pre.bam > NPM1-078.lofreq.bam
/home/programs/lofreq_star-2.1.4_linux-x86-64/bin/lofreq call -b dynamic -C 2 -a 0.00005 -q 20 -Q 20 -m 50 -f /home/reference_genomes/hg19_broad/hg19_all.fasta -l /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed -o NPM1-078.lofreq.vcf NPM1-078.lofreq.bam --no-default-filter
/home/programs/lofreq_star-2.1.4_linux-x86-64/bin/lofreq filter --no-defaults -a 0.0001 -i NPM1-078.lofreq.vcf -o NPM1-078.lofreq.filtered.vcf

# minimap_getitd
/home/programs/samtools-1.7/samtools sort NPM1-078.ErrorCorrectd.bam -o NPM1-078.sorted.bam
/home/programs/samtools-1.7/samtools index NPM1-078.sorted.bam
/home/programs/samtools-1.7/samtools view NPM1-078.sorted.bam -b -h chr13 > NPM1-078.chr13.bam
/usr/bin/bedtools bamtofastq -i NPM1-078.chr13.bam -fq NPM1-078_chr13.fastq
python /home/programs/GET_ITD_1_5_15/getitd//getitd.py -reference /home/programs/GET_ITD_1_5_15/getitd//anno/amplicon.txt -anno /home/programs/GET_ITD_1_5_15/getitd//anno/amplicon_kayser.tsv -forward_adapter AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -reverse_adapter CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -nkern 8 NPM1-078 NPM1-078_chr13.fastq

# coverage_mosdepth
/home/pipelines/Consensus_pipeline_with_espresso/scripts/mosdepth.sh NPM1-078.ErrorCorrectd.bam NPM1-078_umi_cov /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed
/home/pipelines/Consensus_pipeline_with_espresso/scripts/mosdepth.sh NPM1-078.ErrorCorrectd.bam NPM1-078_exon_umi_cov /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/AML_MRD_Exons_sortd.bed

# varscan
/home/programs/samtools-1.7/samtools mpileup -d 1000000 -f /home/reference_genomes/hg19_broad/hg19_all.fasta NPM1-078.ErrorCorrectd.bam > NPM1-078.mpileup
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -jar /home/programs/VarScan.v2.3.9.jar mpileup2snp NPM1-078.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1e-4 --output-vcf 1 > NPM1-078.varscan_snp.vcf
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -jar /home/programs/VarScan.v2.3.9.jar mpileup2indel NPM1-078.mpileup --min-coverage 2 --min-reads2 2 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1e-4 --output-vcf 1 > NPM1-078.varscan_indel.vcf
bgzip -c NPM1-078.varscan_snp.vcf > NPM1-078.varscan_snp.vcf.gz
bgzip -c NPM1-078.varscan_indel.vcf > NPM1-078.varscan_indel.vcf.gz
/home/programs/bcftools-1.9/bcftools index -t NPM1-078.varscan_snp.vcf.gz
/home/programs/bcftools-1.9/bcftools index -t NPM1-078.varscan_indel.vcf.gz
/home/programs/bcftools-1.9/bcftools concat -a NPM1-078.varscan_snp.vcf.gz NPM1-078.varscan_indel.vcf.gz -o NPM1-078.varscan.vcf

# CNS_filegen
/home/programs/samtools-1.7/samtools mpileup -s -x -BQ 0 -q 1 -d 1000000 --skip-indels -f /home/reference_genomes/hg19_broad/hg19_all.fasta NPM1-078.ErrorCorrectd.bam > NPM1-078.mpileup

# These parmeters were taken from the Espresso paper(Abelson et al., Sci. Adv. 2020; 6 : eabe3722)
# A stringent (lower) p value removes all the snps hence, p-value is 1
# /usr/lib/jvm/java-8-openjdk-amd64/bin//java -jar /home/programs/VarScan.v2.3.9.jar pileup2cns NPM1-078".mpileup" --variants SNP --min-coverage 10 --min-reads2 1 --min-avg-qual 30 --min-var-freq 0.0001 --p-value 1 --strand-filter 0 > NPM1-078".cns"

/usr/lib/jvm/java-8-openjdk-amd64/bin//java -jar /home/programs/VarScan.v2.3.9.jar pileup2cns NPM1-078".mpileup" --variants SNP --min-coverage 2 --min-reads2 1 --min-avg-qual 15 --min-var-freq 0.0001 --p-value 1 --strand-filter 0 > NPM1-078".cns"
grep -v '^chrM' NPM1-078".cns" > NPM1-078".nochrM.cns"
/home/pipelines/Consensus_pipeline_with_espresso/scripts/filter_cns.sh NPM1-078".nochrM.cns" NPM1-078".final.cns"

# mutect2_run
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -Xmx10G -jar /home/programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T MuTect2 -R /home/reference_genomes/hg19_broad/hg19_all.fasta -I:tumor NPM1-078.ErrorCorrectd.bam -o NPM1-078.mutect2.vcf --dbsnp /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.vcf -L /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed -nct 25 -contamination 0.02 -mbq 20

# hsmetrics_run
/usr/lib/jvm/java-8-openjdk-amd64/bin//java -jar /home/programs/picard/build/libs/picard.jar CollectHsMetrics I= NPM1-078.ErrorCorrectd.bam O= NPM1-078_hsmetrics.txt BAIT_INTERVALS= /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.interval_list TARGET_INTERVALS= /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.interval_list R= /home/reference_genomes/hg19_broad/hg19_all.fasta VALIDATION_STRINGENCY=LENIENT

# somaticSeq_run
somaticseq_parallel.py --output-directory NPM1-078.somaticseq --genome-reference /home/reference_genomes/hg19_broad/hg19_all.fasta --inclusion-region /home/pipelines/Consensus_pipeline_with_espresso/bedfiles/Probes-XGEN_npm1_cebpa_phf6_sortd.bed --threads 25 --algorithm xgboost  --dbsnp-vcf /home/reference_genomes/dbSNPGATK/dbsnp_138.hg19.somatic.vcf single --bam-file NPM1-078.ErrorCorrectd.bam --mutect2-vcf NPM1-078.mutect2.vcf --vardict-vcf NPM1-078.vardict.vcf --varscan-vcf NPM1-078.varscan.vcf --lofreq-vcf NPM1-078.lofreq.filtered.vcf --sample-name NPM1-078

/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/vcf_sorter.sh NPM1-078.somaticseq/Consensus.sSNV.vcf NPM1-078.somaticseq/somaticseq_snv.vcf
bgzip -c NPM1-078.somaticseq/somaticseq_snv.vcf > NPM1-078.somaticseq/somaticseq_snv.vcf.gz
/home/programs/bcftools-1.9/bcftools index -t NPM1-078.somaticseq/somaticseq_snv.vcf.gz

/home/pipelines/NextSeq_mutation_detector_leukemia/scripts/vcf_sorter.sh NPM1-078.somaticseq/Consensus.sINDEL.vcf NPM1-078.somaticseq/somaticseq_indel.vcf
bgzip -c NPM1-078.somaticseq/somaticseq_indel.vcf > NPM1-078.somaticseq/somaticseq_indel.vcf.gz
/home/programs/bcftools-1.9/bcftools index -t NPM1-078.somaticseq/somaticseq_indel.vcf.gz

/home/programs/bcftools-1.9/bcftools concat -a NPM1-078.somaticseq/somaticseq_snv.vcf.gz NPM1-078.somaticseq/somaticseq_indel.vcf.gz -o NPM1-078.ErrorCorrectd.vcf

perl /home/programs/annovar_latest/annovar//convert2annovar.pl -format vcf4 NPM1-078.ErrorCorrectd.vcf --outfile NPM1-078.ErrorCorrectd.avinput --withzyg --includeinfo

perl /home/programs/annovar_latest/annovar//table_annovar.pl NPM1-078.ErrorCorrectd.avinput --out NPM1-078.somaticseq --remove --protocol refGene,cytoBand,cosmic84,popfreq_all_20150413,avsnp150,intervar_20180118,1000g2015aug_all,clinvar_20170905 --operation g,r,f,f,f,f,f,f --buildver hg19 --nastring '-1' --otherinfo --csvout --thread 10 /home/programs/annovar_latest/annovar//humandb/ --xreffile /home/programs/annovar_latest/annovar//example/gene_fullxref.txt

if [ -s NPM1-078.somaticseq.hg19_multianno.csv ]; then
        python3 /home/pipelines/Consensus_pipeline_with_espresso/scripts/somaticseqoutput-format_v2.py NPM1-078.somaticseq.hg19_multianno.csv NPM1-078.ErrorCorrectd.csv
else
        touch NPM1-078.ErrorCorrectd.csv
fi
