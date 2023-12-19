import pandas as pd
import os, sys
import re

args = sys.argv
sample = args[1]
filepath = args[2]
outfile = args[3]
#coverage = args[4]
#espresso = args[5]
#mutect2 = args[6]

csvfilenames=[filepath+sample+'_final_anno.csv']

writer = pd.ExcelWriter(outfile)
for csvfilename in csvfilenames:
	if os.path.getsize(csvfilename) != 0:
		sheetname=os.path.split(csvfilename)[1]
		df = pd.read_csv(csvfilename)
		print('process file:', csvfilename, 'shape:', df.shape)
		new_sheet_name = os.path.splitext(sheetname)[0]
		new_sheet_name = re.sub (sample,"", new_sheet_name, flags = re.IGNORECASE)
		df.to_excel(writer,sheet_name=new_sheet_name, index=False)
writer.save()
