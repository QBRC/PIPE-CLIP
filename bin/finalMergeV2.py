#!/usr/bin/python
# programmer : bbc
# usage:

import sys
import re
import random
import string
import math
from pysam import *
from pybedtools import BedTool
import rpy2.robjects as robject
from rpy2.robjects.packages import importr
from rpy2.robjects import FloatVector
import argparse as ap

stats = importr('stats')

def fisherTest(p):
	'''First p is cluster p, others are mutation p'''
	R = robject.r
	fps = []
	for i in p[1:]:
		xsq = 2*(p[0]+i)
		fp =  R.pchisq(xsq,**{'df':4,'lower.tail':False,'log.p':True})[0]
		fps.append(-fp)
	return max(fps)

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
	mutationLoc = {}
	kms = {}
	pvalues = {}
	for item in bed_merge:
		if item[7]!=".":
			name = "\t".join(item[0:6])
			if overlap.has_key(name):
				overlap[name].append(item[10])
				mutationLoc[name].append(item[8])
				kms[name].append("("+str(item[14])+","+str(item[15])+")")
				pvalues[name].append(float(item[11]))
			else:
				overlap[name]=[item[10]]
				mutationLoc[name] = [item[8]]
				kms[name]=["("+str(item[14])+","+str(item[15])+")"]
				pvalues[name]=[float(item[6]),float(item[11])]
	print "#chr\tstart\tend\tname\treadsCount\tstrand\tmutation_locations\tmutation_km\t-log(p)"
	for k in overlap.keys():
		fisher_p = fisherTest(pvalues[k])
		#overlap[k].insert(4,str(fisher_p))
		print "%s\t%s\t%s\t%f" % (k,','.join(mutationLoc[k]),','.join(kms[k]),fisher_p)

if __name__=="__main__":
	main()
