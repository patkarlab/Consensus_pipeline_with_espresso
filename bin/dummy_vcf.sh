#!/usr/bin/env bash

lowarm=$1
moletag_count=$2
output_vcf=$3
orientation=$4

orient=$(echo ${orientation} | awk '{print tolower($0)}')
if [[ "${orient}" == "f" ]]
then
	frac=1.1
else
	frac=0.0
fi	


echo "##fileformat=VCFv4.1" > ${output_vcf}
echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> ${output_vcf}
echo "##FORMAT=<ID=ALT,Number=1,Type=Integer,Description=\"Number of alts found\">" >> ${output_vcf}
echo "##FORMAT=<ID=TOT,Number=1,Type=Integer,Description=\"total number of molecular tags at this location\">" >> ${output_vcf}
echo "##FORMAT=<ID=FRAC,Number=1,Type=Float,Description=\"fraction of alt alleles compared to total molecular tags at this location\">" >> ${output_vcf}
echo "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1" >> ${output_vcf}
let "lowarm++"
echo -e "$CHR\t${lowarm}\t.\tX\tX\t.\tPASS\t.\tGT:ALT:TOT:FRAC\t1:0:${moletag_count}:${frac}" >> ${output_vcf}
