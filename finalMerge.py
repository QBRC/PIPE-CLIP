#!/usr/bin/python
# programmer : bbc
# usage:
#add test
import sys
import re
import random
import string
import math
import rpy2.robjects as robject
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector
from pysam import *
from pybedtools import BedTool
import argparse as ap

stats = importr('stats')

def main():
	try:
		cluster_bed = BedTool(sys.argv[1])
	except IOError,message:
		print >> sys.stderr,"Cannot open BAM file.",message
		sys.exit(1)
	
	try:
		mutation_bed = BedTool(sys.argv[2])
	except IOError,message:
		print >> sys.stderr,"Cannot open BAM file.",message
		sys.exit(1)

	bed_merge = cluster_bed.intersect(mutation_bed,wao=True)
	overlap = {}
	for item in bed_merge:
		if item[7]!=".":
			name = item[0]+"\t"+item[1]+"\t"+item[2]+"\t"+item[3]+"\t"+item[4]+"\t"+item[5]+"\t"+item[6]
			if overlap.has_key(name):
				overlap[name].append(item[10])
			else:
				overlap[name]=[item[10]]
	print "#chr\tstart\tent\tname\tread_numer\tstrand\t-log(q)\tmutation_in_cluster"
	for k in overlap.keys():
		print "%s\t%s" % (k,','.join(overlap[k]))

if __name__=="__main__":
	main()
