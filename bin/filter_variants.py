#!/usr/bin/env python3
# This script will take a list of variants and remove them from the input file

import argparse
import pandas as pd

def filter_variants(arguments):
	outfile_name = arguments.output
	input_file_list = arguments.input_file
	variant_list = arguments.variant_list

	df = pd.read_csv(input_file_list)
	df2 = pd.read_csv(variant_list, sep='\t')
	merged_df = (pd.merge(df, df2, left_on=['Chr','Start','Ref','Alt'], right_on=[df2.columns[0], df2.columns[1], df2.columns[2], df2.columns[3]], 
				how='left', indicator=True).query('_merge == "left_only"').drop(columns=['_merge'] + list(df2.columns)))
	merged_df.to_csv(outfile_name, index=False)

def parse_arguments():
	parser = argparse.ArgumentParser(description="Deletion of repetitive variants")		
	parser.add_argument('--output', help='csv output after removal of variants in variant list eg: temp.csv')	# positional argument
	parser.add_argument('--input_file', help='annotated list of variants')
	parser.add_argument('--variant_list', help='list of variants to be filtered')
	return parser.parse_args()

def main ():
	args = parse_arguments()
	filter_variants(args)

if __name__ == "__main__":
	main()