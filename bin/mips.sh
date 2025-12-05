#!/usr/bin/bash

sample=$1
txt_file=$2
sam_orientation=$3
mip_count_file=$4

while read full_line
do
	line=$(echo $full_line | cut -d' ' -f 1);
	ranges=$(echo $full_line | cut -d' ' -f 2);
	lowarm=$(echo $full_line | cut -d' ' -f 3);
	higharm=$(echo $full_line | cut -d' ' -f 4);
	CHR=$(echo $ranges | cut -d':' -f 1);
	NT=$(echo $ranges | cut -d':' -f 2);
	NT=$(echo $NT | cut -d'-' -f 1);
	CHRT=$CHR$'\t'

	mkdir ${line}
	GREPSTR="LC_ALL=C fgrep $CHRT ${sample}_${CHR}${sam_orientation}.sam | grep '$NT' | sort -t: -k 8,8  > ${line}/${line}.$NT.sam"
	eval $GREPSTR

	if [ -s "${line}/${line}.$NT.sam" ]
	then
		perl remove_singles.pl ${line}/${line}.$NT.sam ${mip_count_file}
	fi
	
	perl add_RX_field.pl ${line}/${line}.$NT.sam.nosingles
	cat dict.sam ${line}/${line}.$NT.sam.nosingles.RX > ${line}/${line}.$NT.sam.nosingles.RX_withheader.sam
done < ${txt_file}