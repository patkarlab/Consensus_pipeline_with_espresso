#!/usr/bin/env python
# This script will calculate the probability of a vaf based on the beta matrix

import argparse
import pandas as pd
from scipy.stats import beta
import math

def PVAL(x, alpha_val, beta_val):
	if math.isnan(alpha_val) or math.isnan(beta_val):
		one_minus_pbeta = -1	# Assign -1 for empty alpha / beta values
	else:
		pbeta = beta.cdf(x, alpha_val, beta_val)
		one_minus_pbeta = 1 - pbeta
	return one_minus_pbeta

def BACKGROUND(alpha_value, beta_value):
	if math.isnan(alpha_value) or math.isnan(beta_value):
		background = -1
	else:
		mean = (alpha_value/(alpha_value + beta_value))
		std_dev = mean * (((1-mean) / alpha_value ) ** 0.5)
		background = (mean + (3 * std_dev)) * 100	# Print the % instead of vaf
	return background

def pval_beta(arguments):
	print ("beta matrix file is:", arguments.matrix_file)
	print ("input csv file is:", arguments.input_csv_file)
	print ("output csv file is", arguments.output_csv_file)

	#sheet_to_modify = 'Variants'	
	incsv_file = arguments.input_csv_file
	matrix_file = arguments.matrix_file
	PVAL_column = "PVAL(1-pbeta)"
	background_column = "BACKGROUND"

	df = pd.read_csv(incsv_file)
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
	merged_df = pd.merge(df, df_matrix, left_on=[df.columns[0], df.columns[1]], right_on=[df_matrix.columns[0], df_matrix.columns[1]], how='left')
	merged_df[PVAL_column] = float(-1)		# Introduce an extra column for P value
	merged_df[background_column] = float(-1)	# Another extra column for Background

	for index, row in merged_df.iterrows():
		if row['Ref'] != '-':
			vaf = row['VAF%']
			VAF = (vaf / 100)
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

	merged_df = merged_df[ df.columns.to_list() + [PVAL_column,background_column]]

	# Write back to the SAME file (overwrite)
	merged_df.to_csv(arguments.output_csv_file, index=False)

def parse_arguments():
	parser = argparse.ArgumentParser(description="Calculation of beta distribution matrix")
	parser.add_argument('--matrix_file', help='matrix file with the alpha beta values for each position')	# positional argument
	parser.add_argument('--input_csv_file', help='input csv file with variants')		# 
	parser.add_argument('--output_csv_file', help='output csv file with Error values')	#
	return parser.parse_args()

def main ():
	args = parse_arguments()
	pval_beta(args)

if __name__ == "__main__":
	main()
