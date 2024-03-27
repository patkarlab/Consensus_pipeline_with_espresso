#!/usr/bin/bash 

infile=$1
outfile=$2

grep "^#" ${infile} > ${outfile}
grep -v "^#" ${infile} | awk 'BEGIN{FS=OFS="\t"}$7="PASS"' >> ${outfile}
