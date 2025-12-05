#!/usr/bin/nextflow

process VARCALL {
	tag "${Sample}"
	label 'process_medium'
	input:
		tuple val(Sample), path(mpileup_files), path(cons_sam)
		path(mip_file)
		val(orientation)
	output:
		tuple val(Sample), path("${Sample}_${orientation}_combine_all_vcfs_sorted.vcf")
	script:
	"""
	sam_file_ext=.sam.nosingles.RX_withheader_con.sam
	for files in ${mpileup_files}
	do
		file_name=\$( basename \${files} .mpileup)
		sam_file_name=\${file_name}\${sam_file_ext}
		mip_name=\$(echo \${file_name} | cut -d '.' -f 1)	# Extracting the mip name
		call_variants_from_mpileup.pl \${files} '0' \${file_name}.vcf
		grep chr \${file_name}.vcf | sort -k 2,2n -k 4,4 -k 5,5 > \${file_name}.unmodded
		lines_in_vcf=\$(awk 'END{print NR}' \${file_name}.unmodded )

		if [ -s "\${sam_file_name}" ]
		then
			moletag_count=\$(awk 'END{print NR}' \${sam_file_name})
			moletag_count=\$(echo \${moletag_count} | awk '{print \$1-3}')	#three header lines
		else
			moletag_count=0
		fi
		echo "moletag count \${file_name} for was \$moletag_count"

		# Extracting the lowarm and higharm values
		lowarm=\$( grep -w \$mip_name ${mip_file} | cut -d' ' -f 3)
		higharm=\$( grep -w \$mip_name ${mip_file} | cut -d' ' -f 4)

		if [ \${lines_in_vcf} == "0" ]
		then
			echo "No variants found for \${file_name}.unmodded, \${moletag_count} tagcount"
			#create dummy vcf to keep moletag_count
			dummy_vcf.sh \${lowarm} \${moletag_count} \${file_name}.final.vcf ${orientation}
		else
			simplify_vcfs.pl \${lowarm} \${higharm} \${file_name}.unmodded \${moletag_count} \${file_name}.final.vcf
		fi
	done
	grep -h chr *.final.vcf | sort -k 1,1 -k 2,2n -k 4,4 -k 5,5 > ${Sample}_${orientation}_combine_all_vcfs_sorted.vcf
	"""
}