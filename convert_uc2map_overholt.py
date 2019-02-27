#!/usr/bin/env python
""" 
Written by Will Overholt
6-25-2016
"""

import sys,re,os
from collections import OrderedDict

UCFileName = sys.argv[1]
OTU_DICT = OrderedDict()
count = 0
with open(UCFileName) as UC:
	for line in UC:
		line = line.rstrip()
		Field = line.split("\t")
		if Field[0] == 'S':
			#this is a seed sequence for the dict key
			OTU_DICT.setdefault(Field[8])
		elif Field[0] == 'H':
			if OTU_DICT.has_key(Field[9]):
				try:
					OTU_DICT[Field[9]].append(Field[8])
				except AttributeError:
					OTU_DICT[Field[9]] = [Field[8]]
			#this is a hit against a seed
			#OTU_DICT.setdefault(Field[9]).append(Field[8])
			#OTU_DICT[Field[9]].append(Field[8])

for i,OTU in enumerate(OTU_DICT, 1):
	sys.stdout.write("Derep_OTU{0}\t{1}".format(i, OTU))
	if OTU_DICT[OTU]:
		for s in OTU_DICT[OTU]:
			sys.stdout.write("\t" + s)
		sys.stdout.write("\n")
	else:
		sys.stdout.write("\n")
