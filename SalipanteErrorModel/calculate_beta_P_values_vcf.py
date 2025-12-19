#!/usr/bin/env python
# This script will calculate the probability of a vaf based on the beta matrix

import argparse
import sys
import linecache
import statistics
import pandas as pd
from scipy.stats import beta

def PVAL(x, alpha_val, beta_val):
	pbeta = beta.cdf(x, alpha_val, beta_val)
	one_minus_pbeta = 1 - pbeta
	return one_minus_pbeta

def BACKGROUND(alpha_value, beta_value):
	mean = (alpha_value/(alpha_value + beta_value))
	std_dev = mean * (((1-mean) / alpha_value ) ** 0.5)
	background = mean + (3 * std_dev)
	return background

def pval_beta(arguments):
	print ("beta matrix file is:", arguments.matrix_file)
	print ("input excel file is:", arguments.input_excel_file)
	print ("output excel file is", arguments.output_excel_file)

	sheet_to_modify = 'Variants'	
	excel_file = arguments.input_excel_file
	matrix_file = arguments.matrix_file
	PVAL_column = "PVAL(1-pbeta)"
	background_column = "BACKGROUND"

	df = pd.read_excel(excel_file, sheet_name=None)
	df_matrix = pd.read_csv(matrix_file, sep='\t')

	# Column indices for A/T/G/C
	A_alpha_col = df_matrix.columns[2]
	A_beta_col = df_matrix.columns[3]
	T_alpha_col = df_matrix.columns[4]
	T_beta_col = df_matrix.columns[5]
	G_alpha_col = df_matrix.columns[6]
	G_beta_col = df_matrix.columns[7]
	C_alpha_col = df_matrix.columns[8]
	C_beta_col = df_matrix.columns[9]
	merged_df = pd.merge(df[sheet_to_modify], df_matrix, left_on=[df[sheet_to_modify].columns[0], df[sheet_to_modify].columns[1]], right_on=[df_matrix.columns[0], df_matrix.columns[1]], how='left')
	merged_df[PVAL_column] = int(-1)		# Introduce an extra column for P value
	merged_df[background_column] = int(-1)	# Another extra column for Background

	for index, row in merged_df.iterrows():
		if row['Ref'] != '-':
			vaf_list = [row['VAF%_Mutect2'], row['VAF%_VarScan2'], row['VAF%_VarDict']]
			for vafs in vaf_list:	# Take the first vaf value which is not -1
				if vafs != -1:
					VAF = (vafs / 100)
					break
				else:
					pass
			if row['Alt'] == "A":
				alpha = row[A_alpha_col]
				beta = row[A_beta_col]
				merged_df.loc[index, PVAL_column] = PVAL(VAF, alpha, beta)
				merged_df.loc[index, background_column] = BACKGROUND(alpha, beta)
			elif row['Alt'] == "T":
				alpha = row[T_alpha_col]
				beta = row[T_beta_col]
				merged_df.loc[index, PVAL_column] = PVAL(VAF, alpha, beta)
				merged_df.loc[index, background_column] = BACKGROUND(alpha, beta)
			elif row['Alt'] == "G":
				alpha = row[G_alpha_col]
				beta = row[G_beta_col]
				merged_df.loc[index, PVAL_column] = PVAL(VAF, alpha, beta)
				merged_df.loc[index, background_column] = BACKGROUND(alpha, beta)
			elif row['Alt'] == "C":
				alpha = row[C_alpha_col]
				beta = row[C_beta_col]
				merged_df.loc[index, PVAL_column] = PVAL(VAF, alpha, beta)
				merged_df.loc[index, background_column] = BACKGROUND(alpha, beta)
			else:
				pass
		else:
			pass

	merged_df = merged_df[ df[sheet_to_modify].columns.to_list() + [PVAL_column, background_column]]

	# Write back to the SAME file (overwrite)
	with pd.ExcelWriter(arguments.output_excel_file , engine="openpyxl", mode="w") as writer:
		merged_df.to_excel(writer, sheet_name=sheet_to_modify, index=False)
		for sheet_name, df in df.items():
			if sheet_name != sheet_to_modify:	# Not printing the original sheet
				df.to_excel(writer, sheet_name=sheet_name, index=False)

def parse_arguments():
	parser = argparse.ArgumentParser(description="Calculation of beta distribution matrix")
	parser.add_argument('--matrix_file', help='matrix file with the alpha beta values for each position')	# positional argument
	parser.add_argument('--input_excel_file', help='list of files to make the beta matrix')		# 
	parser.add_argument('--output_excel_file', help='list of files to make the beta matrix')	#
	return parser.parse_args()

def main ():
	args = parse_arguments()
	pval_beta(args)

if __name__ == "__main__":
	main()