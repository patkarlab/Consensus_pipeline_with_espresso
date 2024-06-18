#!/usr/bin/bash 

infile=$1
outfile=$2

grep "^#" ${infile} > ${outfile}
grep -v "^#" ${infile} | sort -k1,1V -k2,2g | awk 'BEGIN{FS=OFS="\t"}$7="PASS"' >> ${outfile}
