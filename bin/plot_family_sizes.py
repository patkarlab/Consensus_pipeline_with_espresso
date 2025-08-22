#!/usr/bin/env python
# This script will take the tag-family-sizes.txt from fgbio 'GroupReadsByUmi' as input

import  matplotlib.pyplot as plt
import sys

family_size_file = sys.argv[1]
output_file = sys.argv[2]

x = []
y = []
with open (family_size_file) as infile:
	header = next(infile)
	for lines in infile:
		family_size = int (lines.strip().split('\t')[0])
		count = float (lines.strip().split('\t')[2]) * 100
		x.append(family_size)
		y.append(count)
		#print (family_size, count)
plt.plot(x,y)
plt.xlim(0, 100)
plt.xlabel('family size')
plt.ylim(0,100)
plt.ylabel('% of total reads')
plt.savefig(f'{output_file}.png', bbox_inches='tight')