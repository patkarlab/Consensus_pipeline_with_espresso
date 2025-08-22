#!/usr/bin/env python3
# This script will combine csv files into one excel file

import pandas as pd
import os, sys

output_excel=sys.argv[1]	# .xlsx file
ErrorCorrCSV=sys.argv[2]
coverage=sys.argv[3]

writer = pd.ExcelWriter(output_excel)

df1 = pd.DataFrame()
mutect2_sheet_name="merged"
if os.path.getsize(ErrorCorrCSV) != 0:
	df1 = pd.read_csv(ErrorCorrCSV)
df1.to_excel(writer,sheet_name=mutect2_sheet_name, index=False)

df2 = pd.DataFrame()
coverage_sheet_name="coverage"
if os.path.getsize(coverage) != 0:
	df2 = pd.read_csv(coverage, sep='\t', header=None)	
df2.to_excel(writer,sheet_name=coverage_sheet_name, index=False, header=False)

writer.save()
