#!/usr/bin/env python3
# This script will combine csv files into one excel file

import pandas as pd
import os, sys

output_excel=sys.argv[1]	# .xlsx file
mutect2=sys.argv[2]
varscan=sys.argv[3]
vardict=sys.argv[4]
coverage=sys.argv[5]

writer = pd.ExcelWriter(output_excel)

df1 = pd.read_csv(mutect2)
mutect2_sheet_name="mutect2"
df1.to_excel(writer,sheet_name=mutect2_sheet_name, index=False)

df2 = pd.read_csv(varscan)
varscan_sheet_name="varscan"
df2.to_excel(writer,sheet_name=varscan_sheet_name, index=False)

df3 = pd.read_csv(vardict)
vardict_sheet_name="vardict"
df3.to_excel(writer,sheet_name=vardict_sheet_name, index=False)

df4 = pd.read_csv(coverage, sep='\t', header=None)
coverage_sheet_name="coverage"
df4.to_excel(writer,sheet_name=coverage_sheet_name, index=False, header=False)

writer.save()
