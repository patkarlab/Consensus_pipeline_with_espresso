#!/usr/bin/env python

import matplotlib.backends.backend_pdf
from matplotlib import pyplot as plt
import csv
import math
import sys

infile = sys.argv[1]	# input file is output from bedtools coverage -d 
outfile = sys.argv[2]	# 

pdf = matplotlib.backends.backend_pdf.PdfPages(outfile)
coverage_val = []
location = []
chrom_start_stop = dict ()
plot_title = dict ()
with open (infile) as input:
	file_handle = csv.reader(input, delimiter="\t")
	for lines in file_handle:
		chrom = lines[0]
		start = int (lines[1])
		end = int(lines[2])
		position = int(lines[-2]) - 1
		coverage = math.log(int(lines[-1]) + 1, 10)
		diff = end - start
		
		location.append(position)
		coverage_val.append(coverage)	
		if int(lines[-2]) == diff:						
			plt.figure (figsize=(6, 5))
			title = lines[-3] + '\n' + str(chrom) + ':' + str(start) + '-' + str(end)
			plt.title(title)
			plt.scatter(location, coverage_val, s=10)
			plt.ylim (0, 5)
			plt.ylabel('(log10)')
			plt.plot (location, coverage_val)
			pdf.savefig()
			plt.close()
			location = []
			coverage_val = []
pdf.close()
