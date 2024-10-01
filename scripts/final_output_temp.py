#!/usr/bin/env python3
# This script will combine CSV files into one Excel file

import pandas as pd
import os
import sys
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Combine CSV files into one Excel file.')
    parser.add_argument('output_excel', type=str, help='Path to the output Excel file (.xlsx)')
    parser.add_argument('ErrorCorrCSV', type=str, help='Path to the error corrected CSV file')
    parser.add_argument('coverage', type=str, help='Path to the coverage file (tab-separated)')
    return parser.parse_args()

args = parse_arguments()

# Assign arguments to variables
output_excel = args.output_excel
ErrorCorrCSV = args.ErrorCorrCSV
coverage = args.coverage

# Create an Excel writer object
writer = pd.ExcelWriter(output_excel)

# Read the first CSV file and write it to an Excel sheet
df1 = pd.read_csv(ErrorCorrCSV)
mutect2_sheet_name = "merged"
df1.to_excel(writer, sheet_name=mutect2_sheet_name, index=False)

# Read the second CSV file (tab-separated) and write it to another Excel sheet
df2 = pd.read_csv(coverage, sep='\t', header=None)
coverage_sheet_name = "coverage"
df2.to_excel(writer, sheet_name=coverage_sheet_name, index=False, header=False)

# Save the Excel file
writer.save()

