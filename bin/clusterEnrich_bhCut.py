#!/usr/bin/python
# programmer : bbc
# usage:


import sys
import re
import random
import string
import pysam
from pysam import *
import argparse as ap
from pybedtools import BedTool
import copy
import clusterEnrich_nb as ce
import math
from collections import Counter


def prepare_argparser():
	description = "Looking for reliable clusters"
	epilog = "For command line options of each command, type %(prog)s COMMAND -h"
	argparser = ap.ArgumentParser(description = description, epilog = epilog)
	argparser.add_argument("-i","--input",dest = "input",type = str, required = True, help = "input file")
	argparser.add_argument("-o","--output",dest = "output",type = str, required = True, help = "output file prefix")
	argparser.add_argument("-f",dest = "fdr",type = float, help = "FDR cutoff", default = 0.05)

	return(argparser)
		


def muEvaluate(mapfile,threshold):
	reliableList = []
	km_p = {}#store km and corresponding p value
	pvalues = []
	total_test = 0
	for row in mapfile:
		fields = row.rstrip().split("\t")
		if(fields[0]!="chrom") and (fields[6] != "NA"):
			pvalues.append(float(fields[6]))
			readLen = str(int(fields[2])-int(fields[1]))
			k = readLen+"_"+fields[4]
			km_p[k]=float(fields[6])
			total_test += 1
	pCount = dict(Counter(pvalues))
	#print pCount
	pRank = ce.freqRank(pCount,True)
	pqDic={}
	for i in pRank.keys():
		try:
			p_rank = pRank[i]
			q_value = ce.BH(i,p_rank,total_test)
			pqDic[i]=q_value
		except:
			print >> sys.stderr,"Cannot find p value in dictionary"
			continue
	mapfile.seek(0)
	for record in mapfile:
		buf = record.rstrip().split("\t")
		if(buf[0]!="chrom"):
			name = str(int(buf[2])-int(buf[1]))+"_"+str(buf[4])
			record_p = km_p[name]
			record_q = pqDic[record_p]
			if record_q >= math.log(threshold)*(-1):
			#buf.append(str(record_q))
				buf[6]=str(record_q)
				reliableList.append(buf)

	return reliableList


def main():
	#for item in muEvaluate_check(sys.argv[1],115065019,0.01):
	#	print item
	
	argparser = prepare_argparser()
	args = argparser.parse_args()
	
	
	try:
		infile = open(args.input,"r")
	except IOError,message:
		print >> sys.stderr, "cannot open coverage file",message
		sys.exit(1)
	filename = args.output
	outputfile = open(filename,"wa")
	print >> outputfile,"#chr\tstart\tend\tname\treadCount\tstrand\t-log(q)"
	for reliable_mu in muEvaluate(infile,args.fdr):
		print >>outputfile,'\t'.join(reliable_mu[0:7])
	

if __name__=="__main__":
	main()
