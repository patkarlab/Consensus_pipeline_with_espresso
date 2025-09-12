#!/usr/bin/env python3
# This script will take input files generated from fill_empty_mips.pl
# The output of the script will be a matrix file for calculate_beta_P_values_vcf.pl

import argparse
import sys
import linecache

def count_beta_dist(arguments):
	# print ("samples are", arguments.samples, type(arguments.samples))
	# print ("output to", arguments.output, type(arguments.output))

	# Check if all samples have the same no. of lines else exit
	line_count = 0
	for index, sample_file in enumerate (arguments.samples):
		with open (sample_file) as infile:
			num_lines = sum(1 for _ in infile)	# Counting the no. of lines in the input file
			if ( index == 0):
				line_count = num_lines
			else:
				if (num_lines != line_count):
					print ("No. of lines in", sample_file, "not same as previous file, Stopping the execution")
					sys.exit () 	# Stop the analysis if inputs are not of same length
				else:
					pass

	# The no. of columns are 9 , their names are : chr, pos, ref, tagcount, A, T, G, C , "null"
	# Iterate over line numbers 
	# mean_percent_error_A = []
	# mean_percent_error_T = []
	# mean_percent_error_G = []
	# mean_percent_error_C = []
	total_tag_count = []
	total_A_count = []
	total_T_count = [] 
	total_G_count = []
	total_C_count = []
	for lines in range(2, line_count+1):	# Ignore the first line with column names
		# The following list are for samples
		tag_count_list = []
		a_count_list = []
		t_count_list = []
		g_count_list = []
		c_count_list = []
		for ind, samp_file in enumerate (arguments.samples):
			line_info = linecache.getline(samp_file, lines).split('\t')
			tag_count = int (line_info[3])
			a_count = int (line_info[4])
			t_count = int (line_info[5])
			g_count = int (line_info[6])
			c_count = int (line_info[7])
			tag_count_list.append(tag_count)
			a_count_list.append(a_count)
			t_count_list.append(t_count)
			g_count_list.append(g_count)
			c_count_list.append(c_count)

		total_tag_count.append(sum(tag_count_list))
		total_A_count.append(sum(a_count_list))
		total_T_count.append(sum(t_count_list))
		total_G_count.append(sum(g_count_list))
		total_C_count.append(sum(c_count_list))

	#print (len(total_tag_count), len(total_A_count))

def parse_arguments():
	parser = argparse.ArgumentParser(description="Calculation of beta distribution matrix")
	parser.add_argument('--output', help='name of the output matrix file eg: beta_matrix.txt')	# positional argument
	parser.add_argument('--samples', nargs='*', help='list of files to make the beta matrix')
	return parser.parse_args()

def main ():
	args = parse_arguments()
	count_beta_dist(args)

if __name__ == "__main__":
	main()