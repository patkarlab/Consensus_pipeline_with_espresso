#!/usr/bin/bash

# to create link for .final.cns files 

folder=$1
samplesheet=$2

for samples in `cat ${samplesheet}`
do
		ln -s ${folder}/${samples}.final.cns ./
done
