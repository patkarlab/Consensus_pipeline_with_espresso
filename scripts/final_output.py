#!/usr/bin/env python3
# This script will combine csv files into one excel file

import pandas as pd
import os, sys

output_excel=sys.argv[1]	# .xlsx file
ErrorCorrCSV=sys.argv[2]
coverage=sys.argv[3]
mutect2=sys.argv[4]

writer = pd.ExcelWriter(output_excel)

df1 = pd.read_csv(ErrorCorrCSV)
mutect2_sheet_name="merged"
df1.to_excel(writer,sheet_name=mutect2_sheet_name, index=False)

df2 = pd.read_csv(coverage, sep='\t', header=None)
coverage_sheet_name="coverage"
df2.to_excel(writer,sheet_name=coverage_sheet_name, index=False, header=False)

df3 = pd.read_csv(mutect2)
mutect2_sheet_name = "mutect2"
df3.to_excel(writer,sheet_name=mutect2_sheet_name, index=False)

writer.save()
