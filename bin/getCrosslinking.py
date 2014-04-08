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
	c_p = p[0]
	m_p = min(p[1:])
	#for i in p[1:]:
	try:
		xsq = -2*math.log(c_p * m_p)
		fp =  R.pchisq(xsq,**{'df':4,'lower.tail':False,'log.p':True})[0]
		fps = -1.0*fp
	except:
		fps = 0.0
	return fps

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
		#print item
		if item[7]!=".":
			#print item
			name = "\t".join(item[0:6])
			if overlap.has_key(name):
				overlap[name].append(item[10])
				mutationLoc[name].append(item[8])
				pvalues[name].append(float(item[11]))
			else:
				overlap[name]=[item[10]]
				mutationLoc[name] = [item[8]]
				pvalues[name]=[float(item[6]),float(item[11])]
	print "#chr\tstart\tstop\tname\treads_count\tstrand\t-log(p)\tmutation_in_cluster\tmutation_locations"
	for k in overlap.keys():
		fisher_p = fisherTest(pvalues[k])
		#overlap[k].insert(4,str(fisher_p))
		print "%s\t%f\t%s\t%s" % (k,fisher_p, ','.join(overlap[k]), ','.join(mutationLoc[k]))

if __name__=="__main__":
	main()
