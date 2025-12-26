#!/usr/bin/env python3

## ./fill_empty_mips.py <mip list>  <multiple variant file>
## input is a non-overlapping list of all the mips and a multiple variant file
## the program will fill in all positions with no variants

import argparse
import re
import linecache

def split_data(line):
	line_split = line.split('\t')
	chr = line_split[0]
	pos = int(line_split[1])
	ref = line_split[2]
	tag = int (line_split[3])
	a_count = int (line_split[4])
	t_count = int (line_split[5])
	g_count = int (line_split[6])
	c_count = int (line_split[7])
	return chr,pos,ref,tag,a_count,t_count,g_count,c_count

def chr_mod(cnr_chr):
	cnr_chr = re.sub ("chr","", cnr_chr, flags = re.IGNORECASE)
	cnr_chr = re.sub ("X","23", cnr_chr, flags = re.IGNORECASE)
	cnr_chr = re.sub ("Y","24", cnr_chr, flags = re.IGNORECASE)
	return int (cnr_chr)

def fill_mips(arguments):
	infile = arguments.input_file
	mips_bed = arguments.bedfile
	outfile = str(infile) + str(".filled")
	default_ref = "X"
	line_no = 1
	print ("infile is", infile)
	print ("bedfile", mips_bed)
	print ("output file is",outfile)

	bed_file = open(mips_bed)
	output = open(outfile, "w")
	first_line = linecache.getline(infile, line_no).strip()
	line_no += 1	# Incrementing the line_no variable

	print (first_line, file=output)
	for lines in bed_file:
		region = lines.split()
		chr = region[0]
		start = int (region[1])
		end = int (region[2])

		variant_data = []
		with open (infile) as input:
			for index, variants in enumerate(input):
				if (index > 0):	# ignore the header
					if chr_mod(variants.split('\t')[0]) == chr_mod(chr):
						variant_data.append(variants)

		idx = 0
		for i in range(start, end + 1):
			#if variant_file_line:
			chr_var, pos_var, ref_var, tag_var, a_var, t_var, g_var, c_var = split_data(variant_data[idx])
			#print (i, pos_var)
			default_ref = "X"
			tag_count = tag_var
			A_count = T_count = G_count = C_count = 0
			for ind, values in enumerate(variant_data):
				chr_var, pos_var, ref_var, tag_var, a_var, t_var, g_var, c_var = split_data(values)
				if (i == pos_var) and chr_mod(chr) == chr_mod(chr_var):
					default_ref = ref_var
					tag_count = tag_var 
					A_count = a_var
					T_count = t_var
					G_count = g_var
					C_count = c_var
					idx = ind
					break
				else:
					pass
			print (chr, i, default_ref, tag_count, A_count, T_count, G_count, C_count, chr, file=output, sep='\t')

	bed_file.close()
	output.close()

def parse_arguments():
	parser = argparse.ArgumentParser(description="script to fill in empty mips positions")		
	parser.add_argument('--input_file', help='file with extension _errors_combined.real_removed')
	parser.add_argument('--bedfile', help='list of non overlapping regions in bed format')
	return parser.parse_args()

def main ():
	args = parse_arguments()
	fill_mips(args)

if __name__ == "__main__":
	main()