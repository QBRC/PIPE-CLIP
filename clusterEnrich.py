#!/usr/bin/python
# programmer : bbc
# usage:

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

def prepare_argparser():
	description = "Find enriched clusters"
	epilog = "For command line options of each command, type %(prog)s COMMAND -h"
	argparser = ap.ArgumentParser(description=description,epilog=epilog)
	argparser.add_argument("-i","--input",dest="infile",type=str,required=True,help="input bam file")
	argparser.add_argument("-n",dest="readCount",help="Use reads in cluster to filter. Default=0",type = int, default = 0)
	argparser.add_argument("-f",dest="fdr",help="Use FDR to filter. Default=0.01. Exclusive to -n",type = float)
	return(argparser)

def BH(pvalue):
	pvalue.sort(key=lambda dup:dup[1],reverse=True)
	m = len(pvalue)
	qvalue = []
	for index in range(1,m+1):
		a = math.log(float(m)/index)
		q = a + float(pvalue[index-1][1])
		qva = -1 * q
		qvalue.append((pvalue[index-1][0],qva))
	return qvalue

def filterByRead(myBed,cutoff):
	count = 1
	for item in myBed:
		if int(item[3])>=cutoff: #the reads number in a cluster
			name = "cluster"+str(count)
			st = [item[0],item[1],item[2],name,item[3],item[4]]
			count +=1
			yield st
	return 

def filterByFDR(myBed,cutoff):
	count = 1
	read = []
	length = []
	R = robject.r
	name = []
	for row in myBed:
		if float(row[3])>1:
			length.append(float(row[2])-float(row[1]))
			read.append(float(row[3]))
			name.append(row)#(row.rstrip())
	cdf = []
	read = FloatVector(read)
	length = FloatVector(length)
	robject.globalenv['read'] = read
	fit = R.glm("read~1",family="poisson",offset=R.log(length))
	intercept =  fit.rx('coefficients')[0][0]
	mu = R.exp(intercept)[0]
	for i in range(len(read)):
		p = R.ppois(read[i],mu*length[i],False,True)[0]
		cdf.append((name[i],p))
	cdf_adjust = BH(cdf)
	for index in cdf_adjust:#range(len(name)):
		if float(index[1]>=math.log(cutoff)*(-1)):
			cluster_name = "cluster_"+str(count)
			count += 1
			st=[index[0][0],index[0][1],index[0][2],cluster_name,index[0][3],index[0][4],str(index[1])]
			yield st
	return

def main():
	argparser = prepare_argparser()
	args = argparser.parse_args()
	try:
		mapped_bam = BedTool(args.infile)
	except IOError,message:
		print >> sys.stderr,"Cannot open BAM file.",message
		sys.exit(1)
	
	mapped_bed = mapped_bam.bam_to_bed()
	bed_merge = mapped_bed.merge(s=True, n=True)
	
	print "#chr\tstart\tend\tcluster_name\treads_in_cluster\tstrand\t-log(q)"
	if args.fdr>0:
		for cl in filterByFDR(bed_merge,args.fdr):
			if int(cl[4])>args.readCount:
				print '\t'.join(cl)


if __name__=="__main__":
	main()
