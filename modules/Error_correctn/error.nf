#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ERROR_CAL {
	// publishDir "$PWD/Final_Output/", mode: 'copy'
	input:
		file (ExcelFiles)
	// output:
	// 	file("hsmetrics_probewise.txt") 
	script:
	"""
	for i in ${ExcelFiles}
	do
		echo \${i} >> file_names.txt
	done
	# Extracting the file names for BNC
	bnc=\$( cat ${params.bnc_excel_files} )
	# Extracting the file names for samples
	samples=\$( grep -vi 'bnc' file_names.txt )

	error_correction.py --samples \${samples} --normals \${bnc}
	for i in `cat ${params.input}`
	do 
		if [ -s \${i}_*_mod.xlsx ]; then
			cp \${i}_*_mod.xlsx ${params.outdir}/\${i}/
		fi
	done
	"""
}

process ERROR_MODEL {
	input:
		tuple val (Sample), file(error_correctd_xl)
	output:
		tuple val (Sample), file()
	script:
	""" 

	"""
}