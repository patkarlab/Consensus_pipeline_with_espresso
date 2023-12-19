#!/usr/bin/bash 

infile=$1
outfile=$2

awk '{s=(NR==1)?"GT":"0/1";$0=$0 "\t" s}1' ${infile} > ${outfile}
